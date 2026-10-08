"""Per-cycle mean quality and low (Q<15) share in FASTQ files: cycles 1-5, 6-10, 146-150, 151."""
import sys
import numpy as np

for arg in sys.argv[1:]:
    name, f = arg.split('=')
    rows = []
    with open(f, 'rb') as fh:
        for i, line in enumerate(fh):
            if i % 4 == 3:
                q = np.frombuffer(line.rstrip(b'\n'), dtype=np.uint8).astype(np.int16) - 33
                if len(q) == 151:
                    rows.append(q)
    a = np.vstack(rows)
    def s(lo, hi):
        x = a[:, lo - 1:hi]
        return '%.2f/%.2f%%' % (x.mean(), 100 * (x < 15).mean())
    print(name, 'n=%d' % len(a), 'meanQ/low%%: 1-5 %s  6-10 %s  50-100 %s  146-150 %s  151 %s' % (s(1, 5), s(6, 10), s(50, 100), s(146, 150), s(151, 151)))
