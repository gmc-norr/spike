"""Statistics for the methods slides, from real data in this folder.

validate's allele-fraction rule is re-implemented from src/validate.rs
(constants ALT_ERROR_RATE 0.001, REF_ERROR_RATE 0.01, AF_TAIL 0.005,
AF_POWER 0.99, MIN_ALT_READS 3) and checked against validate's own verdicts.
"""
import csv, json, math

def pmf(n, p):
    step = math.log(p) - math.log(1 - p); ln = n * math.log(1 - p); out = []
    for k in range(n + 1):
        out.append(math.exp(ln)); ln += (math.log(n - k) if n - k > 0 else 0) - math.log(k + 1) + step
    return out

def floor_(n):
    at_least = 1.0
    for k, pr in enumerate(pmf(n, 0.001)):
        if at_least < 0.005: return max(k, 3)
        at_least -= pr
    return n + 1

def verdict(f, k, n):
    p = f * (1 - 0.01) + (1 - f) * 0.001
    P = pmf(n, p); fl = floor_(n)
    if sum(P[fl:]) < 0.99: return "too shallow", None
    ok = []
    for x in range(n + 1):
        at_most = sum(P[:x + 1]); at_least = sum(P[x:])
        ok.append(x >= fl and at_most >= 0.005 and (f >= 1.0 or at_least >= 0.005))
    acc = [x for x in range(n + 1) if ok[x]]
    return ("PASS" if ok[k] else "FAIL"), (min(acc) / n, max(acc) / n)

rows = list(csv.DictReader(open("snv_af_counts.csv")))
out = []
for r in rows:
    f, k, n = float(r["asked"]), int(r["k"]), int(r["n"])
    v, rng = verdict(f, k, n)
    want = "too shallow" if r["graded"] == "no" else "PASS"
    out.append(dict(asked=f, k=k, n=n, frac=k / n, verdict=v, lo=rng[0] if rng else "", hi=rng[1] if rng else "", agrees=(v == want)))
print("re-implemented rule agrees with validate:", sum(o["agrees"] for o in out), "/", len(out))
with open("snv_accept.csv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)

def wilson(k, n, z=1.96):
    c = (k + z * z / 2) / (n + z * z); h = z * math.sqrt(k * (n - k) / n + z * z / 4) / (n + z * z)
    return c - h, c + h
print("crash real 261/23832", [round(100 * x, 3) for x in wilson(261, 23832)], "spike 1/13016", [round(100 * x, 4) for x in wilson(1, 13016)])

ins = list(csv.DictReader(open("insert.csv")))
a = sorted(int(r["tlen"]) for r in ins if r["kind"] == "real"); b = sorted(int(r["tlen"]) for r in ins if r["kind"] == "spike")
import bisect
D = max(abs(bisect.bisect_right(a, x) / len(a) - bisect.bisect_right(b, x) / len(b)) for x in sorted(set(a + b)))
ne = len(a) * len(b) / (len(a) + len(b))
pks = 2 * sum((-1) ** (j - 1) * math.exp(-2 * j * j * D * D * ne) for j in range(1, 101))
print("KS D", round(D, 4), "n", len(a), len(b), "p approx", pks)
