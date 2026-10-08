#!/bin/bash
set -euo pipefail
# Build the full spiked FASTQ pair from the sample's raw FASTQ pair.
#
# Pipelines such as raredisease start from raw FASTQ. This keeps every raw read
# exactly as it is -- bases, qualities, order and header -- except the
# originals spike removed and did not write back (fastq_removed_reads.txt), and
# adds spike's own reads (named SPIKE_..., from R1.fq.gz/R2.fq.gz) at the end,
# with the header style of the raw file's first record.
#
# Usage: bash fastq.sh RAW_R1 RAW_R2 OUT_R1 OUT_R2 [THREADS]
#
# RAW_R1/RAW_R2 must be the sample's full raw FASTQ (every lane, concatenated:
# `cat` joins gzip files) from the run whose BAM spike was given. A read is
# matched by name: its header's first word without `@` and a trailing /1 or
# /2, the name the aligner gave it in the BAM.
#
# The two mates are read together, record by record, and must hold the same
# names in the same order, as an aligner pairs them. No output may be one of
# the inputs. Each output is written beside its name and moved there only when
# both are complete and checked, so a run that fails changes nothing.
if [ $# -lt 4 ]; then
    echo "Usage: bash fastq.sh RAW_R1 RAW_R2 OUT_R1 OUT_R2 [THREADS]" >&2
    exit 2
fi
RAW=("$1" "$2")
OUT=("$3" "$4")
THREADS="${5:-16}"
DIR="$(cd "$(dirname "$0")" && pwd)"

for f in fastq_removed_reads.txt R1.fq.gz R2.fq.gz; do
    if [ ! -f "$DIR/$f" ]; then
        echo "Error: $DIR/$f not found. Re-run spike." >&2
        exit 1
    fi
done

die() {
    echo "Error: $1" >&2
    exit 1
}

fail() {
    echo "Error: $1" >&2
    echo "RAW_R1/RAW_R2 must be the sample's full raw FASTQ (every lane, concatenated), from the run whose BAM spike was given." >&2
    exit 1
}

# Writing an output over an input destroys that input, so a clash stops the
# run before anything is written.
if [ "${RAW[0]}" -ef "${RAW[1]}" ]; then
    die "${RAW[0]} and ${RAW[1]} are the same file; give the sample's R1 file, then its R2 file."
fi
for o in "${OUT[@]}"; do
    for i in "${RAW[@]}" "$DIR/R1.fq.gz" "$DIR/R2.fq.gz" "$DIR/fastq_removed_reads.txt"; do
        if [ "$o" -ef "$i" ]; then
            die "$o and $i are the same file; an output must not be one of fastq.sh's inputs."
        fi
    done
done
canon() { realpath -m -- "$1" 2> /dev/null || printf '%s\n' "$1"; }
if [ "${OUT[0]}" -ef "${OUT[1]}" ] || [ "$(canon "${OUT[0]}")" = "$(canon "${OUT[1]}")" ]; then
    die "${OUT[0]} and ${OUT[1]} are the same file; OUT_R1 and OUT_R2 must be two files."
fi

if command -v pigz > /dev/null; then
    # The two mates are compressed at the same time.
    ZIP=(pigz -p "$(( THREADS > 1 ? THREADS / 2 : 1 ))")
    UNZIP=(pigz -dc)
else
    ZIP=(gzip)
    UNZIP=(gzip -dc)
fi

TMP="$(mktemp -d)"
PART=()
PIDS=()
# A run that does not finish removes only what it made: its temporary folder,
# its unfinished outputs beside OUT_R1/OUT_R2, and its background jobs.
cleanup() {
    for pid in ${PIDS[@]+"${PIDS[@]}"}; do kill "$pid" 2> /dev/null || true; done
    rm -rf "$TMP"
    if [ ${#PART[@]} -gt 0 ]; then rm -f "${PART[@]}"; fi
}
trap cleanup EXIT
for MATE in 1 2; do
    o="${OUT[$((MATE - 1))]}"
    case "$o" in
        */*) d="${o%/*}"; d="${d:-/}" ;;
        *) d=. ;;
    esac
    PART[$((MATE - 1))]="$(mktemp "$d/.${o##*/}.XXXXXX")"
done

# Both raw mates at once: R1 on the input, R2 through a named pipe. The
# records kept go out unchanged, R1's on the output and R2's to a second
# named pipe; the header style of each mate's first record is kept for
# spike's reads.
PAIR=$(cat <<'AWK'
function stop(why) {
    printf "%s\n", why > ENVIRON["SPIKE_WHY"]
    failed = 1
    exit 1
}
function name_of(header,   name) {
    name = header
    sub(/[ \t].*/, "", name)
    sub(/^@/, "", name)
    sub(/\/[12]$/, "", name)
    return name
}
function style_of(header, mate,   space) {
    space = index(header, " ")
    if (space) return substr(header, space)
    if (header ~ /\/[12]$/) return "/" mate
    return ""
}
BEGIN {
    list = ENVIRON["SPIKE_REMOVED"]
    while ((getline name < list) > 0) {
        if (!(name in removed)) want++
        removed[name] = 0
    }
    close(list)
    raw2 = ENVIRON["SPIKE_RAW2"]
    out2 = ENVIRON["SPIKE_OUT2"]
    # Open R2's output now, so the job reading it never waits.
    printf "" > out2
}
{
    if ((getline line2 < raw2) <= 0)
        stop("RAW_R2 ends after line " (NR - 1) " but RAW_R1 goes on; the mates must have the same records.")
    part = NR % 4
    if (part == 1) {
        record++
        if (substr($0, 1, 1) != "@" || substr(line2, 1, 1) != "@")
            stop("record " record " (line " NR "): a FASTQ header starts with @.")
        name = name_of($0)
        if (name != name_of(line2))
            stop("record " record ": RAW_R1 has " name " and RAW_R2 has " name_of(line2) "; the mates are not in the same order.")
        if (NR == 1) {
            printf "%s", style_of($0, 1) > ENVIRON["SPIKE_STYLE1"]
            close(ENVIRON["SPIKE_STYLE1"])
            printf "%s", style_of(line2, 2) > ENVIRON["SPIKE_STYLE2"]
            close(ENVIRON["SPIKE_STYLE2"])
        }
        skip = (name in removed)
        if (skip) {
            if (removed[name]++ == 0) found++
            dropped++
        }
    } else if (part == 2) {
        length1 = length($0)
        length2 = length(line2)
    } else if (part == 3) {
        if (substr($0, 1, 1) != "+" || substr(line2, 1, 1) != "+")
            stop("record " record " (line " NR "): the third line of a FASTQ record starts with +.")
    } else if (length($0) != length1 || length(line2) != length2) {
        stop("record " record ": a sequence and its quality differ in length.")
    }
    if (!skip) {
        print
        print line2 > out2
    }
}
END {
    if (failed) exit 1
    if ((getline line2 < raw2) > 0)
        stop("RAW_R1 ends after line " NR " but RAW_R2 goes on; the mates must have the same records.")
    if (NR % 4)
        stop("the last record ends after line " (NR % 4) " of its 4.")
    if (found != want)
        stop("found " (found + 0) " of the " (want + 0) " originals listed in fastq_removed_reads.txt.")
    if (dropped != found)
        stop("an original listed in fastq_removed_reads.txt is in the raw pair more than once (" dropped " records for " found " names).")
    print record + 0, found + 0 > ENVIRON["SPIKE_COUNT"]
}
AWK
)

# spike's own reads (SPIKE_...), in the raw header style.
ADD_SPIKE=$(cat <<'AWK'
BEGIN {
    file = ENVIRON["SPIKE_STYLE"]
    if ((getline style < file) <= 0) style = ""
    close(file)
}
NR % 4 == 1 {
    name = substr($0, 2)
    sub(/\/[12]$/, "", name)
    mine = (name ~ /^SPIKE_/)
    if (mine) {
        added++
        print "@" name style
        next
    }
}
mine
END { print added + 0 > ENVIRON["SPIKE_COUNT"] }
AWK
)

echo "R1: ${RAW[0]} -> ${OUT[0]}"
echo "R2: ${RAW[1]} -> ${OUT[1]}"
: > "$TMP/style1"
: > "$TMP/style2"
mkfifo "$TMP/raw2" "$TMP/out2"
"${UNZIP[@]}" "${RAW[1]}" > "$TMP/raw2" &
UNZIP2=$!
"${ZIP[@]}" < "$TMP/out2" > "${PART[1]}" &
ZIP2=$!
PIDS=("$UNZIP2" "$ZIP2")
if ! "${UNZIP[@]}" "${RAW[0]}" \
    | SPIKE_REMOVED="$DIR/fastq_removed_reads.txt" SPIKE_RAW2="$TMP/raw2" SPIKE_OUT2="$TMP/out2" \
      SPIKE_STYLE1="$TMP/style1" SPIKE_STYLE2="$TMP/style2" SPIKE_WHY="$TMP/why" \
      SPIKE_COUNT="$TMP/count" LC_ALL=C awk "$PAIR" \
    | "${ZIP[@]}" > "${PART[0]}"; then
    if [ -s "$TMP/why" ]; then
        fail "$(cat "$TMP/why")"
    fi
    die "reading RAW_R1 or writing OUT_R1 failed; the message is above."
fi
wait "$UNZIP2" || die "reading RAW_R2 failed; the message is above."
wait "$ZIP2" || die "writing OUT_R2 failed; the message is above."
PIDS=()
read -r RECORDS FOUND < "$TMP/count"

for MATE in 1 2; do
    "${UNZIP[@]}" "$DIR/R$MATE.fq.gz" \
        | SPIKE_STYLE="$TMP/style$MATE" SPIKE_COUNT="$TMP/added$MATE" LC_ALL=C awk "$ADD_SPIKE" \
        | "${ZIP[@]}" >> "${PART[$((MATE - 1))]}"
done
ADDED1=$(cat "$TMP/added1")
ADDED2=$(cat "$TMP/added2")
if [ "$ADDED1" -ne "$ADDED2" ]; then
    fail "spike's R1.fq.gz and R2.fq.gz hold $ADDED1 and $ADDED2 of its own reads; they must pair up."
fi

# Both complete and checked: publish them, with the mode a new file gets here.
MODE=$(printf '%o' $(( 0666 & ~0$(umask) )))
for MATE in 1 2; do
    chmod "$MODE" "${PART[$((MATE - 1))]}"
    mv -f -T "${PART[$((MATE - 1))]}" "${OUT[$((MATE - 1))]}"
done
echo "Done: read $RECORDS raw pairs, removed $FOUND original pairs, added $ADDED1 of spike's pairs."
