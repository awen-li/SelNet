#!/usr/bin/env python3
import sys, os, glob, math
import h5py
import numpy as np
from pathlib import Path

def find_dataset(h5, candidates):
    """Return the first dataset that exists (search recursively)."""
    found = []
    def visit(name, obj):
        if isinstance(obj, h5py.Dataset):
            found.append(name)
    h5.visititems(lambda n,o: visit(n,o))
    for cand in candidates:
        for f in found:
            # match by suffix or substring
            if f.endswith(cand) or (cand.lower() in f.lower()):
                return f
    return None

def summary(x):
    x = np.asarray(x).ravel()
    return dict(
        mean=float(np.mean(x)),
        q2_5=float(np.quantile(x, 0.025)),
        q97_5=float(np.quantile(x, 0.975)),
        sd=float(np.std(x, ddof=1))
    )

def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <results_dir>", file=sys.stderr)
        sys.exit(1)
    resdir = Path(sys.argv[1])
    files = sorted(resdir.glob("*.hdf5"))
    if not files:
        print(f"No .hdf5 files in {resdir}", file=sys.stderr)
        sys.exit(2)

    print("run_id,file,a_mean,a_q2.5,a_q97.5,b_mean,b_q2.5,b_q97.5,Ne_mean,Ne_q2.5,Ne_q97.5")
    for fpath in files:
        with h5py.File(fpath, "r") as h5:
            # Heuristic candidates commonly used in spatpg builds
            a_name = find_dataset(h5, ["a", "alpha", "/samples/a", "post/a"])
            b_name = find_dataset(h5, ["b", "beta", "/samples/b", "post/b"])
            ne_name = find_dataset(h5, ["Ne", "ne", "N", "Ne_samples", "/samples/Ne"])

            def get_summ(name):
                if name is None: return dict(mean=math.nan, q2_5=math.nan, q97_5=math.nan, sd=math.nan)
                data = h5[name][:]
                return summary(data)

            aS, bS, neS = get_summ(a_name), get_summ(b_name), get_summ(ne_name)
            run_id = fpath.stem

            print(f"{run_id},{fpath.name},{aS['mean']:.6g},{aS['q2_5']:.6g},{aS['q97_5']:.6g},"
                  f"{bS['mean']:.6g},{bS['q2_5']:.6g},{bS['q97_5']:.6g},"
                  f"{neS['mean']:.6g},{neS['q2_5']:.6g},{neS['q97_5']:.6g}")

if __name__ == "__main__":
    main()
