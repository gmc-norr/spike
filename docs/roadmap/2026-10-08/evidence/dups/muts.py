"""Mutation checks for scripts/duplicates: each mutant must turn a test red."""
import os, shutil, subprocess, sys
D = sys.argv[1]
M = [
 ("measure.py", "read.mapping_quality >= COUNTED_MAPQ", "read.mapping_quality > COUNTED_MAPQ"),
 ("measure.py", "r.query_qualities[q] >= COUNTED_BQ", "r.query_qualities[q] > COUNTED_BQ"),
 ("measure.py", "UNCOUNTED = 0x4 | 0x100 | 0x200 | 0x400 | 0x800", "UNCOUNTED = 0x4 | 0x100 | 0x200 | 0x800"),
 ("measure.py", "    pos0 = pos1 - 1\n    out = []", "    pos0 = pos1\n    out = []"),
 ("measure.py", "if r_spk - r_real <= MARGIN:", "if r_spk - r_real < MARGIN:"),
 ("measure.py", "return \"matters\" if leak >= LEAK", "return \"matters\" if leak > LEAK"),
 ("measure.py", "good >= good_share * n", "good > good_share * n"),
 ("measure.py", "def coverage_ok(bam, chrom, pos1, lo=20, hi=45,", "def coverage_ok(bam, chrom, pos1, lo=20, hi=46,"),
 ("measure.py", "        if r.flag & UNCOUNTED:\n            continue\n        n += 1", "        if r.flag & NOT_PRIMARY:\n            continue\n        n += 1"),
 ("measure.py", "if r.mapping_quality >= 20 and r.is_proper_pair:", "if r.mapping_quality >= 20:"),
 ("measure.py", "if counted(r) and r.reference_start >= start1 - 1", "if r.reference_start >= start1 - 1"),
 ("measure.py", "r.reference_start >= start1 - 1 and", "r.reference_start >= start1 - 2 and"),
 ("measure.py", "    return sum(len(r) for r, _ in sites), sum(t for _, t in sites)", "    return sum(len(r) / t for r, t in sites), len(sites)"),
 ("measure.py", "        elif n in source_dups:", "        elif n in replaced and False:"),
 ("draw_events.py", "if not (spaced(lo, hi) and no_call(lo, hi) and clean(lo, hi)):", "if not (no_call(lo, hi) and clean(lo, hi)):"),
 ("draw_events.py", "if not (spaced(lo, hi) and no_call(lo, hi) and clean(lo, hi)):", "if not (spaced(lo, hi) and clean(lo, hi)):"),
 ("draw_events.py", "if not (spaced(lo, hi) and no_call(lo, hi) and clean(lo, hi)):", "if not (spaced(lo, hi) and no_call(lo, hi)):"),
 ("draw_events.py", 'TRANSITION = {"A": "G", "G": "A", "C": "T", "T": "C"}', 'TRANSITION = {"A": "C", "G": "T", "C": "A", "T": "G"}'),
 ("draw_events.py", "            anchor = e[\"start\"] - 1\n", "            anchor = e[\"start\"]\n"),
 ("draw_events.py", "if all(base_ok(p) for p in range(lo, hi + 1, DEL_STEP if width > 1 else 1)):", "if base_ok(lo):"),
 ("control.py", "if not r.is_proper_pair or r.is_unmapped:", "if r.is_unmapped:"),
 ("control.py", "not (r.flag & (0x4 | 0x100 | 0x800 | 0x400))", "not (r.flag & (0x4 | 0x400))"),
 ("control.py", "len(s[\"rep\"]) == 1 and s[\"dups\"]", "len(s[\"rep\"]) >= 1"),
]
caught = 0
for i, (f, old, new) in enumerate(M, 1):
    p = os.path.join(D, "scripts/duplicates", f)
    src = open(p).read()
    assert src.count(old) == 1, (i, f, old)
    open(p, "w").write(src.replace(old, new))
    shutil.rmtree(os.path.join(D, "scripts/duplicates/__pycache__"), ignore_errors=True)
    r = subprocess.run([sys.executable, "-B", "-m", "pytest", "-q", "-x", "-p", "no:cacheprovider",
                        "scripts/duplicates/test_duplicates.py"], cwd=D, capture_output=True, text=True)
    open(p, "w").write(src)
    red = r.returncode != 0 and " failed" in r.stdout and "error" not in r.stdout.split("\n")[-2].lower()
    last = r.stdout.strip().splitlines()[-1]
    print(f"{i:2d} {f}: {'CAUGHT' if red else 'SURVIVED'}  ({last})", flush=True)
    caught += red
print(f"{caught} of {len(M)} caught", flush=True)
