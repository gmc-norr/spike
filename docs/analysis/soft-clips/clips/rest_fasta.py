"""FASTA of the unexplained ('rest') foreign clips of >= 20 bp, for the relaxed realignment."""
import pickle
D = pickle.load(open("clips.pkl", "rb")); C = D["clips"]; F = pickle.load(open("foreign_cls.pkl", "rb"))
with open("rest_ge20.fa", "w") as f:
    n = 0
    for i, k in F:
        if k == "rest" and C[i]["L"] >= 20: f.write(f">{i}\n{C[i]['clip_ref_orient']}\n"); n += 1
print("rest clips >= 20 bp:", n)
