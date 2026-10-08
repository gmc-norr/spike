import re, sys
txt = open(sys.argv[1]).read()
parts = re.split(r"^=== SAMPLE (\d+) t=([0-9.]+)\n", txt, flags=re.M)[1:]
idle_pat = re.compile(r"syscall \(\)|futex|rayon_core::sleep|wait_until|pthread_cond")
keys = r"(spike::[a-z_]+::[A-Za-z_]+|spike::[a-z_]+ \(|<spike::[^>]*>::[a-z_]+)"
for n, t, s in zip(parts[0::3], parts[1::3], parts[2::3]):
    threads = re.split(r"^Thread \d+ ", s, flags=re.M)[1:]
    for th in threads:
        frames = re.findall(r"^#\d+\s+(?:0x[0-9a-f]+ in )?(.*)$", th, flags=re.M)
        if frames and not idle_pat.search(frames[0]):
            sp = [m.group(0) for f in frames for m in [re.search(keys, f)] if m]
            sp = [x for x in sp if x not in ("spike::main (",)]
            print(n, t, " < ".join(sp[:3]), "|", re.sub(r"<[^<>]*>", "", frames[0])[:50])
            break
    else:
        print(n, t, "IDLE")
