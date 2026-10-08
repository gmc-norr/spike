function stop(why) { printf "%s\n", why > ENVIRON["SPIKE_WHY"]; failed = 1; exit 1 }
function name_of(header,   name) { name = header; sub(/[ \t].*/, "", name); sub(/^@/, "", name); sub(/\/[12]$/, "", name); return name }
BEGIN { FS = "\t"; list = ENVIRON["SPIKE_REMOVED"]; while ((getline name < list) > 0) { if (!(name in removed)) want++; removed[name] = 0 }; close(list); out2 = ENVIRON["SPIKE_OUT2"] }
{
    if (NF != 8) stop("fields")
    if (substr($1,1,1) != "@" || substr($5,1,1) != "@") stop("header")
    if (substr($3,1,1) != "+" || substr($7,1,1) != "+") stop("plus")
    if (length($2) != length($4) || length($6) != length($8)) stop("len")
    name = name_of($1)
    if (name != name_of($5)) stop("order")
    if (name in removed) { if (removed[name]++ == 0) found++; dropped++; next }
    print $1 "\n" $2 "\n" $3 "\n" $4
    print $5 "\n" $6 "\n" $7 "\n" $8 > out2
}
END { if (failed) exit 1; if (found != want) stop("found"); if (dropped != found) stop("dup"); print NR, found + 0 > ENVIRON["SPIKE_COUNT"] }
