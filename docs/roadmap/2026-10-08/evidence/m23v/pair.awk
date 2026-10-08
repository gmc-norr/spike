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
