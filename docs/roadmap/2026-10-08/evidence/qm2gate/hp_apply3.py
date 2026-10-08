"""Throwaway Gate B, step 3 (on top of hp_apply.py and hp_apply2.py): 12 read classes, cut finer at
the poor end (0.25, 0.5, 1, 1.5, 2, 3, 5, 8, 20, 40, 60%), with the key's class field widened to 4 bits."""
def edit(path, pairs):
    s = open(path).read()
    for old, new in pairs:
        assert s.count(old) == 1, (path, old[:80])
        s = s.replace(old, new)
    open(path, "w").write(s)
edit("src/quality.rs", [
("pub const READ_CLASSES: usize = 8;", "pub const READ_CLASSES: usize = 12;"),
("const CLASS_QUANTILES: [f64; READ_CLASSES - 1] = [0.01, 0.03, 0.08, 0.2, 0.4, 0.6, 0.8];",
 "const CLASS_QUANTILES: [f64; READ_CLASSES - 1] = [0.0025, 0.005, 0.01, 0.015, 0.02, 0.03, 0.05, 0.08, 0.2, 0.4, 0.6];"),
("""            tag(0) | hp << 36 | lows << 32 | full << 8 | pos << 4 | flag << 3 | class,
            tag(1) | hp << 36 | lows << 32 | two << 8 | pos << 4 | class,
            tag(2) | hp << 36 | one << 8 | pos << 4 | class,
            tag(3) | hp << 36 | one << 8 | pos << 4,
            tag(4) | pos << 4,""", """            tag(0) | hp << 36 | lows << 32 | full << 12 | pos << 8 | flag << 4 | class,
            tag(1) | hp << 36 | lows << 32 | two << 12 | pos << 8 | class,
            tag(2) | hp << 36 | one << 12 | pos << 8 | class,
            tag(3) | hp << 36 | one << 12 | pos << 8,
            tag(4) | pos << 8,"""),
])
print("applied step 3")
