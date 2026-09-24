#!/usr/bin/env bash
#
# validate_pipeline.sh — Cross-sample DEL spike-and-recover validation
#
# Takes known HG002 DELs from the GIAB T2TQ100 truth set, spikes them into
# a clean 1000 Genomes background BAM, runs an SV caller (Delly), and measures
# sensitivity at multiple VAFs.
#
# Usage:
#   bash scripts/validate_pipeline.sh [--background-bam <path>] [--skip-to <step>]
#
# Every tool is taken from PATH (or from $SAMTOOLS, $DELLY, ... if set), and
# every data path defaults to this repository's layout but can be overridden
# with --giab-dir or with the $REFERENCE / $TRUTH_VCF / $BENCH_BED environment
# variables. Inside a git worktree data/giab_hg38 is a dangling symlink, so pass
# --giab-dir there. EXTRA_TOOL_DIRS below is a last-resort fallback that does
# name one developer's directories; PATH and the environment overrides are the
# portable route.
#
# The run ends non-zero when the pipeline produced no usable result (see
# step 9), so this harness can actually fail. The verdict is attributed to the
# spike-in: step 7b calls the same SVs on the *unspiked* background, and step 8
# fails the run unless the spiked BAM recovers at least --min-gain events more
# than that control does.
#
set -euo pipefail

# ─────────────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

# Directories searched for a tool that is not on PATH. Conda/pixi environments
# are usually not activated when this script runs from cron or from an editor.
# The last entry is machine-specific: it is where delly, truvari, bgzip and
# tabix happen to live on the box this harness was developed on. It is a
# convenience fallback, not a requirement -- anywhere else, put the tools on
# PATH or set $DELLY / $TRUVARI / ... instead.
EXTRA_TOOL_DIRS=(
    "${PROJECT_DIR}/.pixi/envs/default/bin"
    "${HOME}/.pixi/bin"
    "${HOME}/.pixi/envs/default/bin"
    "${HOME}/dev/sv_caller/.pixi/envs/default/bin"
)

# find_tool <name> <override> — the override wins, then PATH, then the
# directories above. Falls back to the bare name so check_tool can report it.
find_tool() {
    local name="$1" override="${2:-}" dir resolved
    if [[ -n "$override" ]]; then
        printf '%s' "$override"
        return 0
    fi
    if resolved="$(command -v "$name" 2>/dev/null)"; then
        printf '%s' "$resolved"
        return 0
    fi
    for dir in "${EXTRA_TOOL_DIRS[@]}"; do
        if [[ -x "${dir}/${name}" ]]; then
            printf '%s' "${dir}/${name}"
            return 0
        fi
    done
    printf '%s' "$name"
}

# Tools
SPIKE="${SPIKE:-${PROJECT_DIR}/target/release/spike}"
SAMTOOLS="$(find_tool samtools "${SAMTOOLS:-}")"
BWAMEM2="$(find_tool bwa-mem2 "${BWAMEM2:-}")"
BCFTOOLS="$(find_tool bcftools "${BCFTOOLS:-}")"
DELLY="$(find_tool delly "${DELLY:-}")"
TRUVARI="$(find_tool truvari "${TRUVARI:-}")"
BGZIP="$(find_tool bgzip "${BGZIP:-}")"
TABIX="$(find_tool tabix "${TABIX:-}")"

# Data paths. GIAB_DIR holds the GIAB/T2TQ100 release; the reference lives in
# reference/ and the HG002 truth files in HG002/.
GIAB_DIR="${GIAB_DIR:-${PROJECT_DIR}/data/giab_hg38}"
REFERENCE="${REFERENCE:-}"
TRUTH_VCF="${TRUTH_VCF:-}"
BENCH_BED="${BENCH_BED:-}"

# Output. Not data/validation: that directory is committed, so the old default
# overwrote tracked fixtures on every run (and left half of them stale).
OUTDIR="${OUTDIR:-${PROJECT_DIR}/validation_run}"

# Parameters
THREADS=8
SEED=42
VAFS=(0.5 0.25 0.1)
MIN_DEL_SIZE=500
MAX_DEL_SIZE=50000
CHROM="chr20"
REGION=""          # optional CHROM:START-END, for a small/fast run
MAX_EVENTS=0       # 0 = no cap
MIN_EVENTS=5       # abort if fewer events survive filtering
MIN_RECALL=""      # optional per-VAF recall floor
MIN_GAIN=1         # spiked TP must beat the background control by this many

# Background sample (1000G NYGC, NA18488, YRI, NovaSeq 2x151, ~30x)
CRAM_URL="https://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239491/NA18488.final.cram"
BG_BAM=""  # set via --background-bam or downloaded

# Control
SKIP_TO=0

# ─────────────────────────────────────────────────────────────────────────────
# Argument parsing
# ─────────────────────────────────────────────────────────────────────────────

while [[ $# -gt 0 ]]; do
    case "$1" in
        --background-bam)
            BG_BAM="$2"; shift 2 ;;
        --skip-to)
            SKIP_TO="$2"; shift 2 ;;
        --outdir)
            OUTDIR="$2"; shift 2 ;;
        --threads)
            THREADS="$2"; shift 2 ;;
        --giab-dir)
            GIAB_DIR="$2"; shift 2 ;;
        --region)
            REGION="$2"; shift 2 ;;
        --vafs)
            read -r -a VAFS <<< "$2"; shift 2 ;;
        --max-events)
            MAX_EVENTS="$2"; shift 2 ;;
        --min-events)
            MIN_EVENTS="$2"; shift 2 ;;
        --min-recall)
            MIN_RECALL="$2"; shift 2 ;;
        --min-gain)
            MIN_GAIN="$2"; shift 2 ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --background-bam <path>  Use this BAM as clean background (skip download)"
            echo "  --skip-to <step>         Skip to step N (1-8)"
            echo "  --outdir <dir>           Output directory [default: <repo>/validation_run;"
            echo "                           a directory whose contents git tracks is refused]"
            echo "  --threads <N>            Thread count [default: 8]"
            echo "  --giab-dir <dir>         GIAB data dir [default: <repo>/data/giab_hg38]"
            echo "  --region <chr:beg-end>   Restrict the whole run to this window"
            echo "                           (slices the background BAM; for a fast run)"
            echo "  --vafs \"<v1 v2 ...>\"     Allele fractions [default: ${VAFS[*]}]"
            echo "  --max-events <N>         Use at most N truth DELs (0 = all)"
            echo "  --min-events <N>         Abort if fewer than N DELs survive [default: 5]"
            echo "  --min-recall <F>         Fail the run if any VAF recalls < F"
            echo "  --min-gain <N>           The highest VAF must recover N more truth"
            echo "                           events than the unspiked background control"
            echo "                           [default: 1; 0 only asks it to match the"
            echo "                           control, which is still a gate]"
            echo "  -h, --help               Show this help"
            echo ""
            echo "Environment overrides: SPIKE, SAMTOOLS, BWAMEM2, BCFTOOLS, DELLY,"
            echo "TRUVARI, BGZIP, TABIX, GIAB_DIR, REFERENCE, TRUTH_VCF, BENCH_BED, OUTDIR."
            exit 0 ;;
        *)
            echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

# Fill in the data paths only now: --giab-dir is parsed above.
REFERENCE="${REFERENCE:-${GIAB_DIR}/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta}"
TRUTH_VCF="${TRUTH_VCF:-${GIAB_DIR}/HG002/GRCh38_HG2-T2TQ100-V1.1_chr20.vcf.gz}"
BENCH_BED="${BENCH_BED:-${GIAB_DIR}/HG002/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed}"

[[ "$MIN_GAIN" =~ ^[0-9]+$ ]] \
    || { echo "ERROR: --min-gain must be a non-negative integer (got '$MIN_GAIN')" >&2; exit 1; }

# An out-of-range --skip-to skipped every step and still reached the verdict,
# which then had nothing to judge. 0 is the default (run everything).
{ [[ "$SKIP_TO" =~ ^[0-9]+$ ]] && (( SKIP_TO <= 8 )); } \
    || { echo "ERROR: --skip-to must be a step number from 1 to 8 (got '$SKIP_TO')" >&2; exit 1; }

# Resolve --outdir once, here, against $PWD like every other path on the
# command line. check_outdir_not_tracked used to resolve a relative --outdir
# against $PROJECT_DIR instead: from scripts/, "--outdir ../data/validation"
# made it ask git about the wrong path and let the run write into tracked
# fixtures, while from outside the repository it refused an unrelated
# directory that merely shared a name. realpath -m does not require the
# directory to exist yet.
OUTDIR="$(realpath -m "$OUTDIR")"

REGION_START=0
REGION_END=0
if [[ -n "$REGION" ]]; then
    if [[ ! "$REGION" =~ ^([^:]+):([0-9]+)-([0-9]+)$ ]]; then
        echo "ERROR: --region must look like chr20:61900000-64200000 (got '$REGION')" >&2
        exit 1
    fi
    CHROM="${BASH_REMATCH[1]}"
    REGION_START="${BASH_REMATCH[2]}"
    REGION_END="${BASH_REMATCH[3]}"
fi

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

log() { echo "[$(date '+%H:%M:%S')] $*"; }

fail() { echo "ERROR: $*" >&2; exit 1; }

# Problems that should not stop the run but must make it end non-zero.
FAILURES=()
note_failure() {
    FAILURES+=("$1")
    echo "FAILURE: $1" >&2
}

# A file's identity for a cache key: path, size, mtime. Missing is a value
# too -- a control computed against a truth VCF that has since been rebuilt
# must not be reused.
file_stamp() {
    if [[ -f "$1" ]]; then
        stat -c '%n:%s:%Y' "$1"
    else
        printf '%s:missing' "$1"
    fi
}

check_tool() {
    local name="$1" path="$2" var="$3"
    if [[ -x "$path" ]] || command -v "$path" >/dev/null 2>&1; then
        return 0
    fi
    echo "ERROR: $name not found (tried: $path)" >&2
    echo "       Put it on PATH or run with ${var}=/path/to/${name}" >&2
    return 1
}

# Count the data (non-header) lines of a VCF.
#
# `grep -c` exits 1 when the count is zero, so a bare $(grep -vc '^#' f) trips
# `set -e` and the old $(grep -vc '^#' f || echo 0) produced the two-line string
# "0\n0" — which made every `[[ $n -lt ... ]]` guard below error out and the run
# continue with zero events. Normalise to a single integer here instead.
count_records() {
    local n=""
    if [[ -f "$1" ]]; then
        n="$(grep -vc '^#' "$1" || true)"
    fi
    n="${n:-0}"
    if [[ ! "$n" =~ ^[0-9]+$ ]]; then
        fail "internal: record count for '$1' is not a number: $(printf '%q' "$n")"
    fi
    printf '%s' "$n"
}

# Read one truvari summary.json as "tp fp fn recall precision f1".
#
# TP is base-side (TP-base): TP + FN is then the number of truth events and
# TP/(TP+FN) is the recall printed next to it. TP-comp counts the *calls* that
# matched, which does not add up with the base-side FN in the same row.
read_truvari_summary() {
    python3 - "$1" <<'PY'
import json, sys
d = json.load(open(sys.argv[1]))
tp = d.get('TP-base', d.get('TP', 0))
fp = d.get('FP', 0)
fn = d.get('FN', 0)
recall = d.get('recall', 0) or 0
precision = d.get('precision', 0) or 0
f1 = d.get('f1', 0) or 0
print(f'{tp} {fp} {fn} {recall:.4f} {precision:.4f} {f1:.4f}')
PY
}

# The highest allele fraction in --vafs: the run's best case, and the one the
# background control and the verdict gate are scored against.
highest_vaf() {
    local best="${VAFS[0]}" vaf
    for vaf in "${VAFS[@]}"; do
        if awk -v a="$vaf" -v b="$best" 'BEGIN{exit !(a>b)}'; then best="$vaf"; fi
    done
    printf '%s' "$best"
}

# The verdict's real question: did the *spike-in* change the result?
#
# The truth DELs are common HG002 variants, so a 1000 Genomes background carries
# several of them already and a caller recovers those with no spike-in at all. A
# plain "TP > 0" gate is therefore satisfied by the background alone. Require
# the spiked run to beat the background control by <gain> events instead.
# Returns 0 (beats it), 1 (does not), 2 (inputs are not counts).
beats_background_control() {
    local tp="$1" control_tp="$2" gain="$3"
    [[ "$tp" =~ ^[0-9]+$ ]] || return 2
    [[ "$control_tp" =~ ^[0-9]+$ ]] || return 2
    [[ "$gain" =~ ^[0-9]+$ ]] || return 2
    (( tp >= control_tp + gain ))
}

# Refuse to write into a directory whose contents git tracks. The default used
# to be data/validation, which is committed, so a plain run silently rewrote
# checked-in fixtures -- and left the ones it did not rewrite stale.
check_outdir_not_tracked() {
    command -v git >/dev/null 2>&1 || return 0
    local repo_root
    repo_root=$(git -C "$PROJECT_DIR" rev-parse --show-toplevel 2>/dev/null) || return 0

    # $OUTDIR is absolute by now. A directory outside the repository cannot
    # hold tracked files, and asking git about it only produces an "outside
    # repository" error -- which, swallowed, used to read as "not tracked".
    [[ "$OUTDIR" == "$repo_root" || "$OUTDIR" == "$repo_root"/* ]] || return 0

    local tracked rc=0
    tracked=$(git -C "$PROJECT_DIR" ls-files -- "$OUTDIR" 2>&1) || rc=$?
    [[ "$rc" -eq 0 ]] \
        || fail "git could not say whether --outdir $OUTDIR is tracked: $tracked"
    [[ -z "$tracked" ]] && return 0
    fail "--outdir $OUTDIR holds files git tracks (e.g. $(echo "$tracked" | head -3 | tr '\n' ' ')).
       This run would overwrite them. Point --outdir at a scratch directory."
}

# ─────────────────────────────────────────────────────────────────────────────
# Step 0: Check prerequisites
# ─────────────────────────────────────────────────────────────────────────────

step0_check_prereqs() {
    log "Step 0: Checking prerequisites..."

    local ok=true
    check_tool "spike"    "$SPIKE"    SPIKE    || ok=false
    check_tool "samtools" "$SAMTOOLS" SAMTOOLS || ok=false
    check_tool "bwa-mem2" "$BWAMEM2"  BWAMEM2  || ok=false
    check_tool "bcftools" "$BCFTOOLS" BCFTOOLS || ok=false
    check_tool "delly"    "$DELLY"    DELLY    || ok=false
    check_tool "truvari"  "$TRUVARI"  TRUVARI  || ok=false
    check_tool "bgzip"    "$BGZIP"    BGZIP    || ok=false
    check_tool "tabix"    "$TABIX"    TABIX    || ok=false

    [[ -f "$REFERENCE" ]]     || { echo "ERROR: Reference not found: $REFERENCE" >&2; ok=false; }
    [[ -f "${REFERENCE}.fai" ]] || { echo "ERROR: Reference index not found: ${REFERENCE}.fai" >&2; ok=false; }
    [[ -f "$TRUTH_VCF" ]]    || { echo "ERROR: Truth VCF not found: $TRUTH_VCF" >&2; ok=false; }
    [[ -f "$BENCH_BED" ]]    || { echo "ERROR: Benchmark BED not found: $BENCH_BED" >&2; ok=false; }

    if [[ "$ok" != "true" ]]; then
        echo "" >&2
        echo "Data paths default to ${GIAB_DIR}; pass --giab-dir <dir> (or set" >&2
        echo "GIAB_DIR / REFERENCE / TRUTH_VCF / BENCH_BED) to point elsewhere." >&2
        fail "Prerequisite check failed. Aborting."
    fi

    check_outdir_not_tracked
    mkdir -p "$OUTDIR"
    log "  All prerequisites OK."
}

# ─────────────────────────────────────────────────────────────────────────────
# Step 1: Prepare background BAM
# ─────────────────────────────────────────────────────────────────────────────

# Slice BG_BAM down to $REGION so the rest of the run stays small. spike, the
# aligner, merge.sh and Delly all then work on the slice, and merge.sh's
# read-name check still holds because spike ran on this very file.
# Delly rejects a BAM whose header names contigs the reference FASTA does not
# have -- which is what a background aligned to the _alt analysis set looks like
# next to a no_alt reference. Catch it here rather than an hour later in step 6.
check_background_contigs() {
    local extra
    extra=$("$SAMTOOLS" view -H "$BG_BAM" \
        | awk -v fai="${REFERENCE}.fai" '
            BEGIN { while ((getline line < fai) > 0) { split(line, f, "\t"); have[f[1]] = 1 } }
            /^@SQ/ {
                for (i = 1; i <= NF; i++)
                    if ($i ~ /^SN:/) {
                        n = substr($i, 4)
                        if (!(n in have) && c < 3) { printf "%s ", n; c++ }
                    }
            }')
    if [[ -n "$extra" ]]; then
        fail "Background BAM names contigs the reference does not have (${extra%% }...).
       Delly will refuse it: use a background BAM aligned to $REFERENCE."
    fi
}

# Find the background BAM a previous run left in ${OUTDIR}/background, for
# --skip-to. Prefer a whole BAM over one of this script's own region slices;
# do not guess when there is more than one candidate.
# Prints the path; returns 1 (none) or 2 (ambiguous).
discover_background_bam() {
    local dir="${OUTDIR}/background" f
    local whole=() sliced=()
    # slice_background_to_region only ever writes <name>.<CHROM>_<beg>_<end>.bam
    # for the chromosome under test; a looser <anything>_<n>_<n>.bam pattern
    # also matched a user BAM called e.g. sample.run_1_2.bam.
    local slice_re="\.${CHROM}_[0-9]+_[0-9]+\.bam$"
    [[ -d "$dir" ]] || return 1
    for f in "$dir"/*.bam; do
        [[ -f "$f" ]] || continue
        if [[ "$f" =~ $slice_re ]]; then
            sliced+=("$f")
        else
            whole+=("$f")
        fi
    done
    if [[ "${#whole[@]}" -eq 1 ]]; then printf '%s' "${whole[0]}"; return 0; fi
    if [[ "${#whole[@]}" -gt 1 ]]; then return 2; fi
    if [[ "${#sliced[@]}" -eq 1 ]]; then printf '%s' "${sliced[0]}"; return 0; fi
    if [[ "${#sliced[@]}" -gt 1 ]]; then return 2; fi
    return 1
}

slice_background_to_region() {
    [[ -n "$REGION" ]] || return 0

    # A resumed run (--skip-to) may have found this run's own slice; slicing it
    # again would only write a second copy under a doubled-up name.
    if [[ "$BG_BAM" == *".${CHROM}_${REGION_START}_${REGION_END}.bam" ]]; then
        log "  Background already sliced to $REGION: $BG_BAM"
        return 0
    fi

    local bg_dir="${OUTDIR}/background"
    mkdir -p "$bg_dir"
    local sliced="${bg_dir}/$(basename "${BG_BAM%.*}").${CHROM}_${REGION_START}_${REGION_END}.bam"

    if [[ ! -f "$sliced" || ! -f "${sliced}.bai" ]]; then
        log "  Slicing background to $REGION ..."
        "$SAMTOOLS" view -b -@ "$THREADS" -T "$REFERENCE" "$BG_BAM" "$REGION" -o "$sliced"
        "$SAMTOOLS" index "$sliced"
    fi

    BG_BAM="$sliced"
    local n
    n=$("$SAMTOOLS" view -c "$BG_BAM")
    [[ "$n" -gt 0 ]] || fail "Region $REGION has no reads in the background BAM"
    log "  Background slice: $BG_BAM ($n reads)"
}

step1_prepare_background() {
    log "Step 1: Preparing background BAM..."

    if [[ -n "$BG_BAM" && -f "$BG_BAM" ]]; then
        log "  Using provided background BAM: $BG_BAM"
        # Verify it has chr20 reads
        local chr20_reads
        chr20_reads=$("$SAMTOOLS" idxstats "$BG_BAM" | awk -v c="$CHROM" '$1==c {print $3}')
        if [[ "${chr20_reads:-0}" -lt 1000 ]]; then
            fail "Background BAM has only ${chr20_reads:-0} reads on $CHROM"
        fi
        log "  Background BAM has $chr20_reads reads on $CHROM"
        check_background_contigs
        slice_background_to_region
        return 0
    fi

    local bg_dir="${OUTDIR}/background"
    mkdir -p "$bg_dir"
    BG_BAM="${bg_dir}/NA18488.chr20.bam"

    if [[ -f "$BG_BAM" && -f "${BG_BAM}.bai" ]]; then
        log "  Background BAM already exists: $BG_BAM (skipping download)"
        check_background_contigs
        slice_background_to_region
        return 0
    fi

    log "  Downloading ${CHROM} from NA18488 (1000G NYGC, NovaSeq 2x151, ~30x)..."
    log "  URL: $CRAM_URL"
    log "  This may take 10-30 minutes depending on network speed."

    # Try remote region extraction first (requires CRAM index at same URL)
    if "$SAMTOOLS" view -b -T "$REFERENCE" -@ "$THREADS" \
            "$CRAM_URL" "$CHROM" 2>/dev/null \
        | "$SAMTOOLS" sort -@ "$THREADS" -o "$BG_BAM" - 2>/dev/null; then
        log "  Remote streaming succeeded."
    else
        log "  Remote streaming failed. Downloading full CRAM first..."
        local full_cram="${bg_dir}/NA18488.full.cram"
        local full_crai="${bg_dir}/NA18488.full.cram.crai"

        # Download CRAM and index
        if [[ ! -f "$full_cram" ]]; then
            curl -L -o "$full_cram" "$CRAM_URL"
            curl -L -o "$full_crai" "${CRAM_URL}.crai"
        fi

        # Extract chr20
        "$SAMTOOLS" view -b -T "$REFERENCE" -@ "$THREADS" \
            "$full_cram" "$CHROM" \
            | "$SAMTOOLS" sort -@ "$THREADS" -o "$BG_BAM" -

        # Clean up full CRAM to save space
        rm -f "$full_cram" "$full_crai"
    fi

    "$SAMTOOLS" index "$BG_BAM"

    # Verify
    local chr20_reads
    chr20_reads=$("$SAMTOOLS" idxstats "$BG_BAM" | awk -v c="$CHROM" '$1==c {print $3}')
    log "  Background BAM: $chr20_reads reads on $CHROM"

    local sample_depth
    sample_depth=$("$SAMTOOLS" depth -r "${CHROM}:10000000-10100000" "$BG_BAM" \
        | awk '{s+=$3; n++} END {if(n>0) printf "%.1f", s/n; else print "0"}')
    log "  Sample depth (${CHROM}:10M-10.1M): ${sample_depth}x"

    check_background_contigs
    slice_background_to_region
}

# ─────────────────────────────────────────────────────────────────────────────
# Step 2: Filter truth VCF to het DELs >= 500bp in benchmark regions
# ─────────────────────────────────────────────────────────────────────────────

step2_filter_truth_vcf() {
    log "Step 2: Filtering truth VCF..."

    local truth_dir="${OUTDIR}/truth"
    mkdir -p "$truth_dir"

    local filtered="${truth_dir}/chr20_dels_filtered.vcf"
    local final_vcf="${truth_dir}/chr20_dels_in_benchmark.vcf"
    local bench_bed_chr20="${truth_dir}/chr20_benchmark.bed"

    # Extract chr20 benchmark regions
    grep "^${CHROM}"$'\t' "$BENCH_BED" > "$bench_bed_chr20" || true
    local n_regions
    n_regions=$(wc -l < "$bench_bed_chr20")
    log "  Benchmark BED: $n_regions regions on $CHROM"
    [[ "$n_regions" -gt 0 ]] || fail "No benchmark regions for $CHROM in $BENCH_BED"

    # Extract VCF header
    zcat "$TRUTH_VCF" | grep '^#' > "$filtered"

    # Filter data lines: SVTYPE=DEL, het, size range, no ALT=*, deduplicate
    zcat "$TRUTH_VCF" \
        | grep -v '^#' \
        | awk -F'\t' -v min="$MIN_DEL_SIZE" -v max="$MAX_DEL_SIZE" \
              -v chrom="$CHROM" -v rstart="$REGION_START" -v rend="$REGION_END" '
        {
            # Must be on the chromosome (or window) under test
            if ($1 != chrom) next

            # Must have SVTYPE=DEL
            if ($8 !~ /SVTYPE=DEL/) next

            # Parse SVLEN (may be positive or negative in T2TQ100).
            # POSIX match() + RSTART/RLENGTH, not the 3-argument match() that
            # only gawk has: under mawk that form is a syntax error and every
            # event is dropped, surfacing only as "Only 0 events after
            # filtering".
            if (!match($8, /SVLEN=-?[0-9]+/)) next
            len = substr($8, RSTART + 6, RLENGTH - 6) + 0
            if (len < 0) len = -len
            if (len < min || len > max) next

            # Skip ALT = "*"
            if ($5 == "*") next

            # Require het genotype (field 10, first sub-field before ":")
            split($10, gt, ":")
            g = gt[1]
            if (g != "0|1" && g != "1|0" && g != "0/1" && g != "1/0") next

            # When --region is given, keep only events wholly inside it
            if (rend > 0 && ($2 < rstart || $2 + len > rend)) next

            print
        }' \
        | sort -k1,1 -k2,2n \
        | awk -F'\t' '!seen[$1":"$2]++' \
        >> "$filtered"

    local n_before
    n_before=$(count_records "$filtered")
    log "  After filtering: $n_before het DELs (${MIN_DEL_SIZE}-${MAX_DEL_SIZE}bp)"

    # Intersect with benchmark BED regions
    # Use bcftools view -T for BED intersection
    # `cat >` rather than `cp`: on this project's ZFS home a cp of a freshly
    # written file stalls for ~75 s (and can produce an all-zero copy).
    cat "$filtered" > "${filtered}.tmp"
    "$BGZIP" -f "${filtered}.tmp"
    "$TABIX" -f -p vcf "${filtered}.tmp.gz"
    "$BCFTOOLS" view -T "$bench_bed_chr20" "${filtered}.tmp.gz" > "${final_vcf}.overlapping"
    rm -f "${filtered}.tmp.gz" "${filtered}.tmp.gz.tbi"

    local n_bench
    n_bench=$(count_records "${final_vcf}.overlapping")
    log "  After benchmark intersection: $n_bench het DELs"

    # Drop events that overlap a kept event. spike rejects overlapping events
    # outright (--allow-overlap only downgrades that to a warning, and says the
    # composition is then approximate), and overlapping truth DELs also make
    # truvari's one-to-one matching ambiguous. Keep the first of each cluster.
    grep '^#' "${final_vcf}.overlapping" > "$final_vcf"
    { grep -v '^#' "${final_vcf}.overlapping" || true; } \
        | awk -F'\t' -v cap="$MAX_EVENTS" '
            BEGIN { kept_chrom=""; kept_end=0; n=0 }
            {
                if (!match($8, /SVLEN=-?[0-9]+/)) next
                len = substr($8, RSTART + 6, RLENGTH - 6) + 0
                if (len < 0) len = -len
                # spike compares [POS, POS+SVLEN) — see validate_event_overlaps.
                if ($1 == kept_chrom && $2 < kept_end) next
                if (cap > 0 && n >= cap) next
                kept_chrom = $1; kept_end = $2 + len; n++
                print
            }' \
        >> "$final_vcf"
    rm -f "${final_vcf}.overlapping"

    local n_after
    n_after=$(count_records "$final_vcf")
    if [[ "$MAX_EVENTS" -gt 0 ]]; then
        log "  After dropping overlaps and capping at $MAX_EVENTS: $n_after het DELs"
    else
        log "  After dropping overlaps: $n_after het DELs"
    fi

    if [[ "$n_after" -lt "$MIN_EVENTS" ]]; then
        fail "Only $n_after events after filtering (need $MIN_EVENTS). Check the truth VCF, --region and --max-events."
    fi

    # bgzip + tabix for spike
    "$BGZIP" -f -c "$final_vcf" > "${final_vcf}.gz"
    "$TABIX" -f -p vcf "${final_vcf}.gz"

    log "  Final truth VCF: ${final_vcf}.gz ($n_after events)"
}

# ─────────────────────────────────────────────────────────────────────────────
# Step 3: Run spike injection at each VAF
# ─────────────────────────────────────────────────────────────────────────────

step3_spike_inject() {
    log "Step 3: Running spike injection..."

    local truth_vcf="${OUTDIR}/truth/chr20_dels_in_benchmark.vcf.gz"
    [[ -f "$truth_vcf" ]] || fail "No filtered truth VCF at $truth_vcf (run step 2 first)"

    for vaf in "${VAFS[@]}"; do
        local spike_out="${OUTDIR}/spike_vaf_${vaf}"

        # Skip only when the directory holds everything step 4 needs. Step 4
        # runs spike's own align.sh and merge.sh, and merge.sh needs
        # replaced_reads.txt; a directory with just the FASTQs and the truth VCF
        # (which is exactly what data/validation ships) used to skip here and
        # then abort in step 4 with "re-run step 3" -- a loop, because step 3
        # skipped again.
        local complete=true f
        for f in R1.fq.gz R2.fq.gz truth.vcf align.sh merge.sh replaced_reads.txt; do
            [[ -f "${spike_out}/${f}" ]] || complete=false
        done
        if [[ "$complete" == "true" ]]; then
            log "  VAF=${vaf}: spike output already complete, skipping."
            continue
        fi
        if [[ -d "$spike_out" ]]; then
            log "  VAF=${vaf}: incomplete spike output, re-running spike."
        fi

        log "  VAF=${vaf}: running spike..."
        # --samtools makes the generated align.sh/merge.sh use the same
        # samtools this script resolved, not whatever happens to be on PATH.
        "$SPIKE" \
            --bam "$BG_BAM" \
            --reference "$REFERENCE" \
            --vcf "$truth_vcf" \
            --allele-fraction "$vaf" \
            --seed "$SEED" \
            -t "$THREADS" \
            --flank 10000 \
            --samtools "$SAMTOOLS" \
            -o "$spike_out" \
            2>&1 | tee "${spike_out}.spike.log"

        # Verify output
        if [[ ! -f "${spike_out}/R1.fq.gz" ]]; then
            fail "spike did not produce ${spike_out}/R1.fq.gz"
        fi

        local n_truth
        n_truth=$(count_records "${spike_out}/truth.vcf")
        log "  VAF=${vaf}: spike produced $n_truth events in truth VCF"
        [[ "$n_truth" -gt 0 ]] || fail "spike wrote an empty truth VCF for VAF=${vaf}"
    done
}

# ─────────────────────────────────────────────────────────────────────────────
# Step 4: Align spiked FASTQ and merge with background BAM
# ─────────────────────────────────────────────────────────────────────────────

step4_align() {
    log "Step 4: Aligning spiked FASTQs and merging with background..."

    # Both scripts are the ones spike itself generated, so this harness
    # exercises the documented workflow instead of a second, divergent copy of
    # it. merge.sh removes the originals by read name (replaced_reads.txt), not
    # by event region, and keeps the background BAM's own @RG SM.
    for vaf in "${VAFS[@]}"; do
        local spike_out="${OUTDIR}/spike_vaf_${vaf}"
        local merged_bam="${spike_out}/merged.bam"

        if [[ -f "$merged_bam" && -f "${merged_bam}.bai" ]]; then
            log "  VAF=${vaf}: merged BAM already exists, skipping."
            continue
        fi

        for s in align.sh merge.sh replaced_reads.txt; do
            [[ -f "${spike_out}/${s}" ]] || fail "${spike_out}/${s} not found. Delete ${spike_out} and re-run from step 3."
        done

        # align.sh invokes the aligner by bare name, so make sure the aligner
        # this script resolved is the one it finds.
        log "  VAF=${vaf}: aligning spiked reads (align.sh)..."
        local aligner_path="$PATH"
        [[ "$BWAMEM2" == */* ]] && aligner_path="$(dirname "$BWAMEM2"):$PATH"
        PATH="$aligner_path" bash "${spike_out}/align.sh" "$REFERENCE" "$THREADS"

        log "  VAF=${vaf}: merging into the background BAM (merge.sh)..."
        bash "${spike_out}/merge.sh" "$BG_BAM" "$REFERENCE" "$THREADS"

        [[ -f "$merged_bam" ]] || fail "merge.sh did not produce $merged_bam"

        # Clean up the spike-only intermediate
        rm -f "${spike_out}/sim.bam" "${spike_out}/sim.bam.bai"

        # Quick stats
        "$SAMTOOLS" flagstat "$merged_bam" > "${spike_out}/flagstat.txt"
        local total_reads
        total_reads=$(head -1 "${spike_out}/flagstat.txt" | awk '{print $1}')
        log "  VAF=${vaf}: $total_reads total reads in merged BAM"
    done
}

# ─────────────────────────────────────────────────────────────────────────────
# Step 5: Run spike validate
# ─────────────────────────────────────────────────────────────────────────────

step5_spike_validate() {
    log "Step 5: Running spike validate..."

    for vaf in "${VAFS[@]}"; do
        local spike_out="${OUTDIR}/spike_vaf_${vaf}"

        log "  VAF=${vaf}: validating..."
        # spike validate exits non-zero when any single check fails, which is
        # expected for a cross-sample spike-in; the harness only insists that it
        # produced a parseable report with at least one passing check.
        "$SPIKE" validate \
            --bam "${spike_out}/merged.bam" \
            --truth "${spike_out}/truth.vcf" \
            --reference "$REFERENCE" \
            --json \
            > "${spike_out}/spike_validate.json" \
            2>"${spike_out}/spike_validate.log" || true

        local summary
        summary=$(python3 -c "
import json
d = json.load(open('${spike_out}/spike_validate.json'))
checks = d.get('checks', [])
passed = sum(1 for c in checks if c.get('pass', False))
print(f'{passed} {len(checks)}')
" 2>/dev/null || echo "")

        if [[ -z "$summary" ]]; then
            note_failure "VAF=${vaf}: spike validate produced no parseable JSON (see ${spike_out}/spike_validate.log)"
            continue
        fi

        local passed total
        read -r passed total <<< "$summary"
        log "  VAF=${vaf}: spike validate: ${passed}/${total} checks passed"
        if [[ "$total" -eq 0 || "$passed" -eq 0 ]]; then
            note_failure "VAF=${vaf}: spike validate passed ${passed}/${total} checks"
        fi
    done
}

# ─────────────────────────────────────────────────────────────────────────────
# Step 6: Call SVs with Delly
# ─────────────────────────────────────────────────────────────────────────────

step6_call_svs() {
    log "Step 6: Calling SVs with Delly..."

    for vaf in "${VAFS[@]}"; do
        local spike_out="${OUTDIR}/spike_vaf_${vaf}"
        local delly_bcf="${spike_out}/delly.bcf"
        local delly_vcf="${spike_out}/delly.vcf.gz"

        if [[ -f "$delly_vcf" && -f "${delly_vcf}.tbi" ]]; then
            log "  VAF=${vaf}: Delly output already exists, skipping."
            continue
        fi

        log "  VAF=${vaf}: running delly call..."
        "$DELLY" call \
            -t DEL \
            -g "$REFERENCE" \
            -o "$delly_bcf" \
            "${spike_out}/merged.bam" \
            2>"${spike_out}/delly.log"

        # Convert BCF to VCF.gz
        "$BCFTOOLS" view "$delly_bcf" \
            | "$BGZIP" -c > "$delly_vcf"
        "$TABIX" -f -p vcf "$delly_vcf"

        local n_calls
        n_calls=$("$BCFTOOLS" view -H "$delly_vcf" | wc -l)
        log "  VAF=${vaf}: Delly found $n_calls DEL calls"
    done
}

# ─────────────────────────────────────────────────────────────────────────────
# Step 7: Benchmark with Truvari
# ─────────────────────────────────────────────────────────────────────────────

step7_truvari_bench() {
    log "Step 7: Benchmarking with Truvari..."

    for vaf in "${VAFS[@]}"; do
        local spike_out="${OUTDIR}/spike_vaf_${vaf}"
        local truvari_out="${spike_out}/truvari"
        local truth_bgz="${spike_out}/truth.vcf.gz"

        # Prepare spike truth VCF for truvari (needs contig headers + bgzip + tabix)
        if [[ ! -f "$truth_bgz" ]]; then
            # Add contig headers from reference .fai (truvari requires them)
            local truth_fixed="${spike_out}/truth_fixed.vcf"
            awk '/^##fileformat/' "${spike_out}/truth.vcf" > "$truth_fixed"
            awk -F'\t' '{printf "##contig=<ID=%s,length=%s>\n", $1, $2}' "${REFERENCE}.fai" >> "$truth_fixed"
            grep '^##' "${spike_out}/truth.vcf" \
                | grep -v '^##fileformat' \
                | { grep -v '^##contig=' || true; } >> "$truth_fixed"
            grep '^#CHROM' "${spike_out}/truth.vcf" >> "$truth_fixed"
            grep -v '^#' "${spike_out}/truth.vcf" >> "$truth_fixed"
            "$BGZIP" -c "$truth_fixed" > "$truth_bgz"
            "$TABIX" -f -p vcf "$truth_bgz"
            rm -f "$truth_fixed"
        fi

        # Remove previous truvari output (it requires empty dir)
        rm -rf "$truvari_out"

        log "  VAF=${vaf}: running truvari bench..."
        "$TRUVARI" bench \
            -b "$truth_bgz" \
            -c "${spike_out}/delly.vcf.gz" \
            -o "$truvari_out" \
            -f "$REFERENCE" \
            --passonly \
            -r 500 \
            -p 0.5 \
            -P 0.5 \
            -s "$MIN_DEL_SIZE" \
            2>"${spike_out}/truvari.log" || true

        # Show results
        if [[ -f "${truvari_out}/summary.json" ]]; then
            local tp fp fn recall precision f1
            read -r tp fp fn recall precision f1 \
                < <(read_truvari_summary "${truvari_out}/summary.json")
            log "  VAF=${vaf}: Recall=${recall}  Precision=${precision}  F1=${f1}  TP=${tp}  FP=${fp}  FN=${fn}"
        else
            note_failure "VAF=${vaf}: truvari produced no summary (check ${spike_out}/truvari.log)"
        fi
    done
}

# ─────────────────────────────────────────────────────────────────────────────
# Step 7b: Background control — what the background recovers with no spike-in
# ─────────────────────────────────────────────────────────────────────────────

# Run exactly the same caller and the same benchmark on the *unspiked*
# background BAM. Everything it recovers is a truth DEL the background sample
# carries itself, and no gate built on the spiked numbers alone can tell those
# apart from the spike-in's own contribution.
step7b_background_control() {
    log "Step 7b: Background control (same calls on the unspiked background)..."

    local ctl_dir="${OUTDIR}/background_control"
    mkdir -p "$ctl_dir"

    local best_vaf truth_bgz
    best_vaf="$(highest_vaf)"
    truth_bgz="${OUTDIR}/spike_vaf_${best_vaf}/truth.vcf.gz"
    [[ -f "$truth_bgz" ]] \
        || fail "No spike truth VCF at $truth_bgz to score the control against (run steps 3-7 first)"

    # Step 7 re-runs truvari for every VAF on every run, so a cached control
    # would put a fresh spiked number next to a floor computed from a different
    # background, region or truth set. Key the cache on all three and redo the
    # control whenever any of them moved.
    local key_file="${ctl_dir}/cache_key.txt" key
    key="bg=$(file_stamp "$BG_BAM") region=${REGION:-all} ref=$(file_stamp "$REFERENCE")"
    key="$key truth=$(file_stamp "$truth_bgz") minsize=${MIN_DEL_SIZE}"
    if [[ ! -f "$key_file" || "$(cat "$key_file")" != "$key" ]]; then
        [[ -f "$key_file" ]] && log "  Control inputs changed; recomputing it."
        rm -rf "${ctl_dir}/delly.bcf" "${ctl_dir}/delly.bcf.csi" \
               "${ctl_dir}/delly.vcf.gz" "${ctl_dir}/delly.vcf.gz.tbi" \
               "${ctl_dir}/truvari" "$key_file"
    fi

    if [[ -f "${ctl_dir}/delly.vcf.gz" && -f "${ctl_dir}/delly.vcf.gz.tbi" ]]; then
        log "  Control Delly output already exists, skipping."
    else
        log "  Running delly call on the background alone ($BG_BAM)..."
        "$DELLY" call \
            -t DEL \
            -g "$REFERENCE" \
            -o "${ctl_dir}/delly.bcf" \
            "$BG_BAM" \
            2>"${ctl_dir}/delly.log"
        "$BCFTOOLS" view "${ctl_dir}/delly.bcf" \
            | "$BGZIP" -c > "${ctl_dir}/delly.vcf.gz"
        "$TABIX" -f -p vcf "${ctl_dir}/delly.vcf.gz"
    fi

    if [[ ! -f "${ctl_dir}/truvari/summary.json" ]]; then
        rm -rf "${ctl_dir}/truvari"
        log "  Benchmarking the control against the same truth VCF..."
        "$TRUVARI" bench \
            -b "$truth_bgz" \
            -c "${ctl_dir}/delly.vcf.gz" \
            -o "${ctl_dir}/truvari" \
            -f "$REFERENCE" \
            --passonly \
            -r 500 \
            -p 0.5 \
            -P 0.5 \
            -s "$MIN_DEL_SIZE" \
            2>"${ctl_dir}/truvari.log" || true
    fi

    if [[ -f "${ctl_dir}/truvari/summary.json" ]]; then
        local tp fp fn recall precision f1 n_truth
        read -r tp fp fn recall precision f1 \
            < <(read_truvari_summary "${ctl_dir}/truvari/summary.json")
        n_truth=$(count_records "${OUTDIR}/spike_vaf_${best_vaf}/truth.vcf")
        printf '%s\n' "$key" > "$key_file"
        log "  Background alone recovers ${tp}/${n_truth} truth DELs (recall ${recall}) — this is the floor the spike-in has to beat."
    else
        note_failure "background control produced no truvari summary (check ${ctl_dir}/truvari.log); nothing in this run can be attributed to the spike-in"
    fi
}

# ─────────────────────────────────────────────────────────────────────────────
# Step 8: Summary report
# ─────────────────────────────────────────────────────────────────────────────

step8_summarize() {
    log "Step 8: Generating summary..."

    local summary_tsv="${OUTDIR}/validation_summary.tsv"
    local best_vaf
    best_vaf="$(highest_vaf)"

    printf "VAF\tN_truth\tN_delly\tTP\tFP\tFN\tRecall\tPrecision\tF1\tSpike_validate\n" \
        > "$summary_tsv"

    # The background control goes in first: every row below it has to be read
    # against this floor, not against zero.
    local control_json="${OUTDIR}/background_control/truvari/summary.json"
    local control_tp=""
    if [[ -f "$control_json" ]]; then
        local c_tp c_fp c_fn c_recall c_precision c_f1 c_calls=0 c_truth=0
        read -r c_tp c_fp c_fn c_recall c_precision c_f1 \
            < <(read_truvari_summary "$control_json")
        control_tp="$c_tp"
        c_truth=$(count_records "${OUTDIR}/spike_vaf_${best_vaf}/truth.vcf")
        if [[ -f "${OUTDIR}/background_control/delly.vcf.gz" ]]; then
            c_calls=$("$BCFTOOLS" view -H "${OUTDIR}/background_control/delly.vcf.gz" 2>/dev/null | wc -l)
        fi
        printf "background\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$c_truth" "$c_calls" "$c_tp" "$c_fp" "$c_fn" \
            "$c_recall" "$c_precision" "$c_f1" "n/a" \
            >> "$summary_tsv"
    fi

    for vaf in "${VAFS[@]}"; do
        local spike_out="${OUTDIR}/spike_vaf_${vaf}"
        local truvari_out="${spike_out}/truvari"

        local n_truth=0 n_calls=0
        local tp=0 fp=0 fn=0 recall="N/A" precision="N/A" f1="N/A"
        local spike_pass="N/A"

        # Truth count
        n_truth=$(count_records "${spike_out}/truth.vcf")

        # Delly call count
        if [[ -f "${spike_out}/delly.vcf.gz" ]]; then
            n_calls=$("$BCFTOOLS" view -H "${spike_out}/delly.vcf.gz" 2>/dev/null | wc -l)
        fi

        # Truvari results
        if [[ -f "${truvari_out}/summary.json" ]]; then
            read -r tp fp fn recall precision f1 \
                < <(read_truvari_summary "${truvari_out}/summary.json")
            # The verdict gate. Recovering events the background already
            # carries validates nothing, so the highest VAF has to beat the
            # background control by --min-gain events. (A bare "TP > 0" gate
            # passes on the background alone: on chr20:61.9-64.2 Mb the
            # unspiked NA18488 background already yields TP=3 of 8.)
            if [[ "$vaf" == "$best_vaf" ]]; then
                local rc=0
                beats_background_control "$tp" "${control_tp:-}" "$MIN_GAIN" || rc=$?
                if [[ "$rc" -eq 2 ]]; then
                    note_failure "VAF=${vaf} (highest): no usable background control (control TP='${control_tp:-<none>}', TP='$tp'), so nothing here can be attributed to the spike-in"
                elif [[ "$rc" -ne 0 ]]; then
                    note_failure "VAF=${vaf} (highest): truvari recovered $tp of $n_truth truth DELs but the unspiked background alone recovers $control_tp — the spike-in contributed $((tp - control_tp)), need at least $MIN_GAIN (--min-gain)"
                fi
            fi
            if [[ -n "$MIN_RECALL" ]] \
                && awk -v r="$recall" -v m="$MIN_RECALL" 'BEGIN{exit !(r<m)}'; then
                note_failure "VAF=${vaf}: recall $recall below --min-recall $MIN_RECALL"
            fi
        elif [[ "$vaf" == "$best_vaf" ]]; then
            # The gate above sits inside the branch, so a highest VAF with no
            # truvari output at all used to skip it: the row printed N/A and
            # the run still said VALIDATION PASSED. Step 7's own complaint
            # about a missing summary only covers the case where step 7 ran,
            # and --skip-to 8 lands here having measured nothing.
            note_failure "VAF=${vaf} (highest): no truvari summary at ${truvari_out}/summary.json, so this run measured nothing that could be attributed to the spike-in"
        fi

        # Spike validate results
        if [[ -f "${spike_out}/spike_validate.json" ]]; then
            spike_pass=$(python3 -c "
import json
try:
    d = json.load(open('${spike_out}/spike_validate.json'))
    checks = d.get('checks', [])
    passed = sum(1 for c in checks if c.get('pass', False))
    print(f'{passed}/{len(checks)}')
except Exception:
    print('N/A')
" 2>/dev/null || echo "N/A")
        fi

        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$vaf" "$n_truth" "$n_calls" "$tp" "$fp" "$fn" \
            "$recall" "$precision" "$f1" "$spike_pass" \
            >> "$summary_tsv"
    done

    echo ""
    echo "╔══════════════════════════════════════════════════════════════════╗"
    echo "║            SPIKE VALIDATION PIPELINE — RESULTS                 ║"
    echo "╚══════════════════════════════════════════════════════════════════╝"
    echo ""
    column -t -s $'\t' "$summary_tsv"
    echo ""
    echo "Full results:  ${OUTDIR}"
    echo "Summary TSV:   ${summary_tsv}"
    echo ""
}

# ─────────────────────────────────────────────────────────────────────────────
# Step 9: Verdict
# ─────────────────────────────────────────────────────────────────────────────

step9_verdict() {
    if [[ "${#FAILURES[@]}" -eq 0 ]]; then
        log "VALIDATION PASSED"
        return 0
    fi

    echo "" >&2
    echo "VALIDATION FAILED — ${#FAILURES[@]} problem(s):" >&2
    local f
    for f in "${FAILURES[@]}"; do
        echo "  - $f" >&2
    done
    exit 1
}

# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

main() {
    log "=== Spike Validation Pipeline ==="
    log "Output: $OUTDIR"
    log "VAFs:   ${VAFS[*]}"
    log "Chrom:  $CHROM"
    [[ -n "$REGION" ]] && log "Region: $REGION"
    log "DEL size: ${MIN_DEL_SIZE}-${MAX_DEL_SIZE}bp"
    echo ""

    step0_check_prereqs

    if [[ "$SKIP_TO" -le 1 ]]; then
        step1_prepare_background
    else
        # Skipped step 1: discover an existing background BAM.
        if [[ -z "$BG_BAM" || ! -f "$BG_BAM" ]]; then
            local found_bg rc=0
            found_bg="$(discover_background_bam)" || rc=$?
            case "$rc" in
                0) ;;
                2) fail "Several BAMs in ${OUTDIR}/background; pass --background-bam to say which one." ;;
                *) fail "No background BAM in ${OUTDIR}/background. Run without --skip-to first, or pass --background-bam." ;;
            esac
            BG_BAM="$found_bg"
            log "  Discovered background BAM: $BG_BAM"
        fi
        check_background_contigs
        slice_background_to_region
    fi

    [[ "$SKIP_TO" -le 2 ]] && step2_filter_truth_vcf
    [[ "$SKIP_TO" -le 3 ]] && step3_spike_inject
    [[ "$SKIP_TO" -le 4 ]] && step4_align
    [[ "$SKIP_TO" -le 5 ]] && step5_spike_validate
    [[ "$SKIP_TO" -le 6 ]] && step6_call_svs
    [[ "$SKIP_TO" -le 7 ]] && step7_truvari_bench
    # Always: the verdict is meaningless without the floor it is measured from.
    step7b_background_control
    step8_summarize

    log "=== Pipeline complete ==="
    step9_verdict
}

# Running the file executes the pipeline; sourcing it only defines the helpers,
# which is how the unit tests in src/main.rs exercise them.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    main "$@"
fi
