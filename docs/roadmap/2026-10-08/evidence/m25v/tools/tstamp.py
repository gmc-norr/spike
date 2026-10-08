import subprocess, sys, time, os
out = sys.argv[1]
cmd = sys.argv[2:]
t0 = time.monotonic()
p = subprocess.Popen(cmd, stderr=subprocess.PIPE, stdout=subprocess.DEVNULL, env=dict(os.environ, RUST_LOG="info"))
with open(out, "w") as f:
    for line in p.stderr:
        f.write("%.3f %s" % (time.monotonic() - t0, line.decode(errors="replace")))
    rc = p.wait()
    ru = os.wait4 if False else None
    f.write("%.3f EXIT %d\n" % (time.monotonic() - t0, rc))
