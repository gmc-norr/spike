function stop(why) { printf "%s\n", why > ENVIRON["WHY"]; failed = 1; exit 1 }
function name_of(header,   name) { name = header; sub(/[ \t].*/, "", name); sub(/^@/, "", name); sub(/\/[12]$/, "", name); return name }
BEGIN { list = ENVIRON["SPIKE_REMOVED"]; while ((getline name < list) > 0) { if (!(name in removed)) want++; removed[name] = 0 }; close(list); names = ENVIRON["NAMES"] }
{
    part = NR % 4
    if (part == 1) {
        record++
        if (substr($0, 1, 1) != "@") stop("header")
        name = name_of($0)
        print name > names
        skip = (name in removed)
        if (skip) { if (removed[name]++ == 0) found++; dropped++ }
    } else if (part == 2) { l = length($0) }
    else if (part == 3) { if (substr($0, 1, 1) != "+") stop("plus") }
    else if (length($0) != l) stop("len")
    if (!skip) print
}
END { if (failed) exit 1; if (NR % 4) stop("trunc"); if (found != want) stop("found"); if (dropped != found) stop("dup"); print record + 0, found + 0 > ENVIRON["COUNT"] }
