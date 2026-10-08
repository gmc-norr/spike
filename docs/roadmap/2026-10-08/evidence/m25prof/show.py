import re, sys
txt = open(sys.argv[1]).read()
want = sys.argv[2]
parts = re.split(r"^(=== SAMPLE \d+ t=[0-9.]+)\n", txt, flags=re.M)[1:]
idle_pat = re.compile(r"syscall \(\)|futex|rayon_core::sleep|wait_until|pthread_cond")
for hdr, s in zip(parts[0::2], parts[1::2]):
    threads = re.split(r"^Thread \d+ ", s, flags=re.M)[1:]
    for t in threads:
        frames = re.findall(r"^#\d+\s+(?:0x[0-9a-f]+ in )?(.*)$", t, flags=re.M)
        if frames and not idle_pat.search(frames[0]):
            j = "\n".join(frames)
            if re.search(want, j):
                print(hdr)
                for fr in frames[:int(sys.argv[3])]:
                    fr = re.sub(r"<[^<>]*>", "<>", fr); fr = re.sub(r"<[^<>]*>", "<>", fr)
                    print("    ", fr[:150])
            break
