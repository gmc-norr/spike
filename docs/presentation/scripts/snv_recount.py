import subprocess, re, csv, math, json
events = [l.strip() for l in open("snv_events.txt")]
val = json.load(open("snv.validate.json"))
obs = [c["observed"] for c in val["checks"] if c["check"] == "allele_freq"]
def base_at(pos, rpos, cigar, seq):
    # pos, rpos 1-based; return base at reference pos or None
    r, q = rpos, 0
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        n = int(n)
        if op in "M=X":
            if r <= pos < r + n:
                return seq[q + pos - r]
            r += n; q += n
        elif op in "IS": q += n
        elif op in "DN":
            if r <= pos < r + n: return None
            r += n
    return None
rows = []
for e, o in zip(events, obs):
    _, chrom, pos, ref, rest = e.split(":")
    alt, af = rest.split(";af=")
    pos = int(pos)
    out = subprocess.run(["samtools", "view", "-q", "20", "-F", "0xF04", "snv/merged.bam", f"{chrom}:{pos}-{pos}"],
                         capture_output=True, text=True).stdout
    votes = {}
    for line in out.splitlines():
        f = line.split("\t")
        b = base_at(pos, int(f[3]), f[5], f[9])
        if b is None or b.upper() not in "ACGT": continue
        votes.setdefault(f[0], set()).add(b.upper())
    counts = {"A":0,"C":0,"G":0,"T":0}
    for name, s in votes.items():
        if len(s) == 1: counts[next(iter(s))] += 1
    n = sum(counts.values()); k = counts[alt]
    frac = k / n if n else float("nan")
    z = 1.96
    if n:
        c = (k + z*z/2) / (n + z*z); h = z*math.sqrt(k*(n-k)/n + z*z/4) / (n + z*z)
        lo, hi = c - h, c + h
    else:
        lo = hi = float("nan")
    m = re.match(r"too shallow \(([\d.]+) at (\d+) reads", o)
    vfrac = float(m.group(1)) if m else float(o)
    vn = int(m.group(2)) if m else None
    rows.append(dict(asked=float(af), k=k, n=n, frac=round(frac, 4), lo=round(lo, 4), hi=round(hi, 4),
                     validate=vfrac, validate_n=vn, graded="no" if m else "yes", match=abs(round(frac, 2) - vfrac) < 0.006 and (vn is None or vn == n)))
for r in rows: print(r)
print("all match validate:", all(r["match"] for r in rows), sum(r["match"] for r in rows), "/", len(rows))
with open("snv_af_counts.csv", "w") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
