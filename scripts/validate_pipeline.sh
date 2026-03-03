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
set -euo pipefail

# ─────────────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
GIAB_DIR="${PROJECT_DIR}/data/giab_hg38"

# Tools
SPIKE="${PROJECT_DIR}/target/release/spike"
SAMTOOLS="/home/parlar_ai/.pixi/bin/samtools"
BWAMEM2="/home/parlar_ai/.pixi/bin/bwa-mem2"
BCFTOOLS="/home/parlar_ai/dev/sv_caller/.pixi/envs/default/bin/bcftools"
DELLY="/home/parlar_ai/dev/sv_caller/.pixi/envs/default/bin/delly"
TRUVARI="/home/parlar_ai/dev/sv_caller/.pixi/envs/default/bin/truvari"
BGZIP="/home/parlar_ai/dev/sv_caller/.pixi/envs/default/bin/bgzip"
TABIX="/home/parlar_ai/dev/sv_caller/.pixi/envs/default/bin/tabix"

# Data paths
REFERENCE="${GIAB_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
TRUTH_VCF="${GIAB_DIR}/GRCh38_HG2-T2TQ100-V1.1_chr20.vcf.gz"
BENCH_BED="${GIAB_DIR}/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed"

# Output
OUTDIR="${PROJECT_DIR}/data/validation"

# Parameters
THREADS=8
SEED=42
VAFS=(0.5 0.25 0.1)
MIN_DEL_SIZE=500
MAX_DEL_SIZE=50000
CHROM="chr20"

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
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --background-bam <path>  Use this BAM as clean background (skip download)"
            echo "  --skip-to <step>         Skip to step N (1-8)"
            echo "  --outdir <dir>           Output directory [default: data/validation]"
            echo "  --threads <N>            Thread count [default: 8]"
            echo "  -h, --help               Show this help"
            exit 0 ;;
        *)
            echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

log() { echo "[$(date '+%H:%M:%S')] $*"; }

check_tool() {
    local name="$1" path="$2"
    if [[ ! -x "$path" ]]; then
        echo "ERROR: $name not found at $path" >&2
        return 1
    fi
}

# ─────────────────────────────────────────────────────────────────────────────
# Step 0: Check prerequisites
# ─────────────────────────────────────────────────────────────────────────────

step0_check_prereqs() {
    log "Step 0: Checking prerequisites..."

    local ok=true
    check_tool "spike"    "$SPIKE"    || ok=false
    check_tool "samtools" "$SAMTOOLS" || ok=false
    check_tool "bwa-mem2" "$BWAMEM2"  || ok=false
    check_tool "bcftools" "$BCFTOOLS" || ok=false
    check_tool "delly"    "$DELLY"    || ok=false
    check_tool "truvari"  "$TRUVARI"  || ok=false
    check_tool "bgzip"    "$BGZIP"    || ok=false
    check_tool "tabix"    "$TABIX"    || ok=false

    [[ -f "$REFERENCE" ]]     || { echo "ERROR: Reference not found: $REFERENCE" >&2; ok=false; }
    [[ -f "${REFERENCE}.fai" ]] || { echo "ERROR: Reference index not found: ${REFERENCE}.fai" >&2; ok=false; }
    [[ -f "$TRUTH_VCF" ]]    || { echo "ERROR: Truth VCF not found: $TRUTH_VCF" >&2; ok=false; }
    [[ -f "$BENCH_BED" ]]    || { echo "ERROR: Benchmark BED not found: $BENCH_BED" >&2; ok=false; }

    if [[ "$ok" != "true" ]]; then
        echo "Prerequisite check failed. Aborting." >&2
        exit 1
    fi

    mkdir -p "$OUTDIR"
    log "  All prerequisites OK."
}

# ─────────────────────────────────────────────────────────────────────────────
# Step 1: Prepare background BAM
# ─────────────────────────────────────────────────────────────────────────────

step1_prepare_background() {
    log "Step 1: Preparing background BAM..."

    if [[ -n "$BG_BAM" && -f "$BG_BAM" ]]; then
        log "  Using provided background BAM: $BG_BAM"
        # Verify it has chr20 reads
        local chr20_reads
        chr20_reads=$("$SAMTOOLS" idxstats "$BG_BAM" | awk -v c="$CHROM" '$1==c {print $3}')
        if [[ "$chr20_reads" -lt 1000 ]]; then
            echo "ERROR: Background BAM has only $chr20_reads reads on $CHROM" >&2
            exit 1
        fi
        log "  Background BAM has $chr20_reads reads on $CHROM"
        return 0
    fi

    local bg_dir="${OUTDIR}/background"
    mkdir -p "$bg_dir"
    BG_BAM="${bg_dir}/NA18488.chr20.bam"

    if [[ -f "$BG_BAM" && -f "${BG_BAM}.bai" ]]; then
        log "  Background BAM already exists: $BG_BAM (skipping download)"
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
    grep "^${CHROM}	" "$BENCH_BED" > "$bench_bed_chr20" || true
    local n_regions
    n_regions=$(wc -l < "$bench_bed_chr20")
    log "  Benchmark BED: $n_regions regions on $CHROM"

    # Extract VCF header
    zcat "$TRUTH_VCF" | grep '^#' > "$filtered"

    # Filter data lines: SVTYPE=DEL, het, size range, no ALT=*, deduplicate
    zcat "$TRUTH_VCF" \
        | grep -v '^#' \
        | awk -F'\t' -v min="$MIN_DEL_SIZE" -v max="$MAX_DEL_SIZE" '
        {
            # Must have SVTYPE=DEL
            if ($8 !~ /SVTYPE=DEL/) next

            # Parse SVLEN (may be positive or negative in T2TQ100)
            match($8, /SVLEN=(-?[0-9]+)/, a)
            len = a[1] + 0
            if (len < 0) len = -len
            if (len < min || len > max) next

            # Skip ALT = "*"
            if ($5 == "*") next

            # Require het genotype (field 10, first sub-field before ":")
            split($10, gt, ":")
            g = gt[1]
            if (g != "0|1" && g != "1|0" && g != "0/1" && g != "1/0") next

            print
        }' \
        | sort -k1,1 -k2,2n \
        | awk -F'\t' '!seen[$1":"$2]++' \
        >> "$filtered"

    local n_before
    n_before=$(grep -vc '^#' "$filtered" || echo 0)
    log "  After filtering: $n_before het DELs (${MIN_DEL_SIZE}-${MAX_DEL_SIZE}bp)"

    # Intersect with benchmark BED regions
    # Use bcftools view -T for BED intersection
    cp "$filtered" "${filtered}.tmp"
    "$BGZIP" -f "${filtered}.tmp"
    "$TABIX" -f -p vcf "${filtered}.tmp.gz"
    "$BCFTOOLS" view -T "$bench_bed_chr20" "${filtered}.tmp.gz" > "$final_vcf"
    rm -f "${filtered}.tmp.gz" "${filtered}.tmp.gz.tbi"

    local n_after
    n_after=$(grep -vc '^#' "$final_vcf" || echo 0)
    log "  After benchmark intersection: $n_after het DELs"

    if [[ "$n_after" -lt 5 ]]; then
        echo "ERROR: Only $n_after events after filtering. Check truth VCF format." >&2
        exit 1
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

    for vaf in "${VAFS[@]}"; do
        local spike_out="${OUTDIR}/spike_vaf_${vaf}"

        if [[ -f "${spike_out}/R1.fq.gz" && -f "${spike_out}/truth.vcf" ]]; then
            log "  VAF=${vaf}: spike output already exists, skipping."
            continue
        fi

        log "  VAF=${vaf}: running spike..."
        "$SPIKE" \
            --bam "$BG_BAM" \
            --reference "$REFERENCE" \
            --vcf "$truth_vcf" \
            --allele-fraction "$vaf" \
            --seed "$SEED" \
            -t "$THREADS" \
            --flank 10000 \
            -o "$spike_out" \
            2>&1 | tee "${spike_out}.spike.log"

        # Verify output
        if [[ ! -f "${spike_out}/R1.fq.gz" ]]; then
            echo "ERROR: spike did not produce ${spike_out}/R1.fq.gz" >&2
            exit 1
        fi

        local n_truth
        n_truth=$(grep -vc '^#' "${spike_out}/truth.vcf" || echo 0)
        log "  VAF=${vaf}: spike produced $n_truth events in truth VCF"
    done
}

# ─────────────────────────────────────────────────────────────────────────────
# Step 4: Align spiked FASTQ and merge with background BAM
# ─────────────────────────────────────────────────────────────────────────────

step4_align() {
    log "Step 4: Aligning spiked FASTQs and merging with background..."

    # Build a BED of all event regions (±flank) for excluding from background
    local events_bed="${OUTDIR}/truth/event_regions.bed"
    # Use the spike truth VCF from the first VAF run (has END/SVLEN in spike format)
    local first_truth="${OUTDIR}/spike_vaf_${VAFS[0]}/truth.vcf"
    if [[ -f "$first_truth" ]]; then
        grep -v '^#' "$first_truth" \
            | awk -F'\t' -v flank=10000 '{
                # Parse END from INFO
                match($8, /END=([0-9]+)/, a)
                end = a[1] + 0
                if (end == 0) {
                    match($8, /SVLEN=-?([0-9]+)/, b)
                    end = $2 + b[1]
                }
                s = $2 - flank; if (s < 0) s = 0
                e = end + flank
                printf "%s\t%d\t%d\n", $1, s, e
            }' \
            | sort -k1,1 -k2,2n \
            | awk -F'\t' '
                BEGIN {chr=""; s=0; e=0}
                {
                    if ($1 != chr || $2 > e) {
                        if (chr != "") printf "%s\t%d\t%d\n", chr, s, e
                        chr=$1; s=$2; e=$3
                    } else { if ($3 > e) e=$3 }
                }
                END { if (chr != "") printf "%s\t%d\t%d\n", chr, s, e }
            ' > "$events_bed"
    else
        echo "ERROR: No spike truth VCF found at $first_truth" >&2
        exit 1
    fi

    local n_regions
    n_regions=$(wc -l < "$events_bed")
    log "  Event regions BED: $n_regions merged intervals"

    for vaf in "${VAFS[@]}"; do
        local spike_out="${OUTDIR}/spike_vaf_${vaf}"
        local spike_bam="${spike_out}/spike_only.bam"
        local merged_bam="${spike_out}/sim.bam"

        if [[ -f "$merged_bam" && -f "${merged_bam}.bai" ]]; then
            log "  VAF=${vaf}: merged BAM already exists, skipping."
            continue
        fi

        # Step 4a: Align spiked FASTQ
        log "  VAF=${vaf}: aligning spiked reads..."
        "$BWAMEM2" mem -t "$THREADS" \
            -R "@RG\\tID:spike_${vaf}\\tSM:SPIKE\\tPL:ILLUMINA\\tLB:lib1" \
            "$REFERENCE" \
            "${spike_out}/R1.fq.gz" "${spike_out}/R2.fq.gz" \
            2>"${spike_out}/align.log" \
            | "$SAMTOOLS" sort -@ "$THREADS" -o "$spike_bam" -
        "$SAMTOOLS" index "$spike_bam"

        # Step 4b: Extract background reads OUTSIDE event regions
        local bg_outside="${spike_out}/bg_outside.bam"
        log "  VAF=${vaf}: extracting background reads outside event regions..."
        "$SAMTOOLS" view -b -@ "$THREADS" \
            -T "$REFERENCE" \
            -L "$events_bed" -U "$bg_outside" \
            "$BG_BAM" "$CHROM" > /dev/null
        "$SAMTOOLS" index "$bg_outside"

        # Step 4c: Merge background (outside events) + spiked (event regions)
        log "  VAF=${vaf}: merging into full-chromosome BAM..."
        "$SAMTOOLS" merge -f -@ "$THREADS" "$merged_bam" "$bg_outside" "$spike_bam"
        "$SAMTOOLS" index "$merged_bam"

        # Clean up intermediates
        rm -f "$spike_bam" "${spike_bam}.bai" "$bg_outside" "${bg_outside}.bai"

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
        "$SPIKE" validate \
            --bam "${spike_out}/sim.bam" \
            --truth "${spike_out}/truth.vcf" \
            --reference "$REFERENCE" \
            --json \
            > "${spike_out}/spike_validate.json" \
            2>"${spike_out}/spike_validate.log" || true

        # Show summary
        if [[ -f "${spike_out}/spike_validate.json" ]]; then
            local summary
            summary=$(python3 -c "
import json, sys
try:
    d = json.load(open('${spike_out}/spike_validate.json'))
    checks = d.get('checks', [])
    passed = sum(1 for c in checks if c.get('pass', False))
    print(f'{passed}/{len(checks)} checks passed')
except:
    print('could not parse JSON')
" 2>/dev/null || echo "N/A")
            log "  VAF=${vaf}: spike validate: $summary"
        else
            log "  VAF=${vaf}: spike validate did not produce JSON (check ${spike_out}/spike_validate.log)"
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
            "${spike_out}/sim.bam" \
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

    local bench_bed="${OUTDIR}/truth/chr20_benchmark.bed"

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
            grep '^##' "${spike_out}/truth.vcf" | grep -v '^##fileformat' >> "$truth_fixed"
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
            python3 -c "
import json
d = json.load(open('${truvari_out}/summary.json'))
recall = d.get('recall', 0) or 0
precision = d.get('precision', 0) or 0
f1 = d.get('f1', 0) or 0
tp = d.get('TP-comp', d.get('TP', 0))
fp = d.get('FP', 0)
fn = d.get('FN', 0)
print(f'  Recall={recall:.3f}  Precision={precision:.3f}  F1={f1:.3f}  TP={tp}  FP={fp}  FN={fn}')
"
        else
            log "  VAF=${vaf}: truvari did not produce summary (check ${spike_out}/truvari.log)"
        fi
    done
}

# ─────────────────────────────────────────────────────────────────────────────
# Step 8: Summary report
# ─────────────────────────────────────────────────────────────────────────────

step8_summarize() {
    log "Step 8: Generating summary..."

    local summary_tsv="${OUTDIR}/validation_summary.tsv"

    printf "VAF\tN_truth\tN_delly\tTP\tFP\tFN\tRecall\tPrecision\tF1\tSpike_validate\n" \
        > "$summary_tsv"

    for vaf in "${VAFS[@]}"; do
        local spike_out="${OUTDIR}/spike_vaf_${vaf}"
        local truvari_out="${spike_out}/truvari"

        local n_truth=0 n_calls=0
        local tp=0 fp=0 fn=0 recall="N/A" precision="N/A" f1="N/A"
        local spike_pass="N/A"

        # Truth count
        if [[ -f "${spike_out}/truth.vcf" ]]; then
            n_truth=$(grep -vc '^#' "${spike_out}/truth.vcf" || echo 0)
        fi

        # Delly call count
        if [[ -f "${spike_out}/delly.vcf.gz" ]]; then
            n_calls=$("$BCFTOOLS" view -H "${spike_out}/delly.vcf.gz" 2>/dev/null | wc -l)
        fi

        # Truvari results
        if [[ -f "${truvari_out}/summary.json" ]]; then
            read -r tp fp fn recall precision f1 < <(python3 -c "
import json
d = json.load(open('${truvari_out}/summary.json'))
tp = d.get('TP-comp', d.get('TP', 0))
fp = d.get('FP', 0)
fn = d.get('FN', 0)
recall = d.get('recall', 0) or 0
precision = d.get('precision', 0) or 0
f1 = d.get('f1', 0) or 0
print(f'{tp} {fp} {fn} {recall:.4f} {precision:.4f} {f1:.4f}')
")
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
except:
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
# Main
# ─────────────────────────────────────────────────────────────────────────────

main() {
    log "=== Spike Validation Pipeline ==="
    log "Output: $OUTDIR"
    log "VAFs:   ${VAFS[*]}"
    log "Chrom:  $CHROM"
    log "DEL size: ${MIN_DEL_SIZE}-${MAX_DEL_SIZE}bp"
    echo ""

    step0_check_prereqs

    [[ "$SKIP_TO" -le 1 ]] && step1_prepare_background
    # If we skipped step 1, discover existing background BAM
    if [[ -z "$BG_BAM" || ! -f "$BG_BAM" ]]; then
        local found_bg="${OUTDIR}/background/NA18488.chr20.bam"
        if [[ -f "$found_bg" ]]; then
            BG_BAM="$found_bg"
            log "  Discovered background BAM: $BG_BAM"
        else
            echo "ERROR: No background BAM found. Run without --skip-to first." >&2
            exit 1
        fi
    fi

    [[ "$SKIP_TO" -le 2 ]] && step2_filter_truth_vcf
    [[ "$SKIP_TO" -le 3 ]] && step3_spike_inject
    [[ "$SKIP_TO" -le 4 ]] && step4_align
    [[ "$SKIP_TO" -le 5 ]] && step5_spike_validate
    [[ "$SKIP_TO" -le 6 ]] && step6_call_svs
    [[ "$SKIP_TO" -le 7 ]] && step7_truvari_bench
    step8_summarize

    log "=== Pipeline complete ==="
}

main "$@"
