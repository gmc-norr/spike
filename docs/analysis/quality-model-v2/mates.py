"""Mate crash link from paired FASTQ: P(both crash) / (P(R1) P(R2)), and crash shares. Usage: mates.py PREFIX..."""
import sys
def quals(path):
    with open(path) as f:
        for i, l in enumerate(f):
            if i % 4 == 3:
                yield l.rstrip("\n")
crash = lambda q: len(q) >= 20 and sum(ord(c) - 33 < 15 for c in q[-20:]) >= 10
for pre in sys.argv[1:]:
    n = c1 = c2 = both = 0
    for a, b in zip(quals(pre + "_R1.fq"), quals(pre + "_R2.fq")):
        x, y = crash(a), crash(b); n += 1; c1 += x; c2 += y; both += x and y
    p1, p2, p12 = c1 / n, c2 / n, both / n
    print(f"{pre}: pairs {n}  P(R1) {100*p1:.2f}%  P(R2) {100*p2:.2f}%  P(both) {100*p12:.3f}%  ratio {p12/(p1*p2):.2f}")
