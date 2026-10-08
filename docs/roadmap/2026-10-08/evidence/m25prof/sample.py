import gdb, os, signal, threading, time
out = open(os.environ["PROF_OUT"], "w")
gdb.execute("set pagination off")
gdb.execute("set confirm off")
gdb.execute("handle SIGINT stop print nopass")
gdb.execute("handle SIGPIPE nostop noprint pass")
gdb.execute("starti")
pid = gdb.selected_inferior().pid
t0 = time.time()
def killer():
    while True:
        time.sleep(0.1)
        try:
            os.kill(pid, signal.SIGINT)
        except ProcessLookupError:
            return
threading.Thread(target=killer, daemon=True).start()
n = 0
while True:
    try:
        gdb.execute("continue", to_string=True)
    except gdb.error:
        break
    try:
        bt = gdb.execute("thread apply all bt 40", to_string=True)
    except gdb.error:
        break
    n += 1
    out.write("=== SAMPLE %d t=%.3f\n" % (n, time.time() - t0))
    out.write(bt)
    out.flush()
out.close()
