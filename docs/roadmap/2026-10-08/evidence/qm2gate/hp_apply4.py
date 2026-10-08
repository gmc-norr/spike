"""Throwaway Gate B, step 4 (on top of steps 1-3): SPIKE_QX bit 64 keeps the lows-in-last-16 bin in
the quality context's coarser levels 2-4 too (not only 0-1)."""
def edit(path, pairs):
    s = open(path).read()
    for old, new in pairs:
        assert s.count(old) == 1, (path, old[:80])
        s = s.replace(old, new)
    open(path, "w").write(s)
edit("src/quality.rs", [
("""            tag(2) | hp << 36 | one << 12 | pos << 8 | class,
            tag(3) | hp << 36 | one << 12 | pos << 8,
            tag(4) | pos << 8,""", """            tag(2) | hp << 36 | low_all << 32 | one << 12 | pos << 8 | class,
            tag(3) | hp << 36 | low_all << 32 | one << 12 | pos << 8,
            tag(4) | low_all << 32 | pos << 8,"""),
("""        let tag = |level: u64| level << 58;
        [
            tag(0) | hp << 36""", """        let low_all = if qx() & 64 != 0 { lows_bin(state.recent_low) as u64 } else { 0 };
        let tag = |level: u64| level << 58;
        [
            tag(0) | hp << 36"""),
])
print("applied step 4")
