#!/usr/bin/env python3
"""Parameter-sweep driver for the `ltgs` Luttinger-Tisza solver.

Accepts the same arguments as the `ltgs` binary, but every numeric argument
takes `nargs='+'` with the following meaning:

  * 1 value   ->  held fixed at that value.
  * 3 values  ->  a `seq`-style axis "START STEP STOP" (STOP inclusive).
  * 2 values  ->  error (ambiguous).
  * anything else -> error.

Every swept argument becomes an independent axis; the driver runs `ltgs` over
the full Cartesian product of the axes, accumulates the results, and writes a
single self-describing HDF5 file for a downstream plotting script.

Example (phase-diagram over J2 x Jc1 at fixed J1, L):

    ./ltgs_sweep.py 6 --J1 1 --J2 -1 0.1 1 --Jc1 0 0.25 1 -o sweep.h5

Output datasets (grid_shape = tuple of swept-axis lengths, in `swept_params`
order; trailing axes as noted):

    axes/<name>          1-D values of each swept parameter
    E_per_spin           grid_shape           Luttinger-Tisza energy / spin
    lambda_min           grid_shape           min eigenvalue of M(k)  (= 2*E/spin)
    k_star_vec           grid_shape + (3,)     ordering wavevector k* (rad/coord-unit)
    k_star_idx           grid_shape + (3,)     k* grid index
    sublattice_weights   grid_shape + (n_sl,)  |s_a|^2 of the k* eigenvector

Root attributes: `swept_params` (comma-joined, axis order), `grid_shape`,
`n_sl`, `fixed_<name>` for each fixed parameter, and `ltgs_bin`.
"""
import argparse
import glob
import math
import os
import subprocess
import sys
import tempfile

import h5py
import numpy as np

# (name passed to ltgs, is_integer). Order here fixes the output axis order.
# L is the positional argument; the rest are --flags understood by ltgs
# (declare_LJ123 + declare_DM in include/cli_bits.hpp).
SCALAR_PARAMS = [
    ("L", True),
    ("J1", False),
    ("J2", False),
    ("J3", False),
    ("Jd", False),
    ("Jc1", False),
    ("D", False),
    ("d", False),
]


def expand(name, vals, is_int):
    """Turn a raw nargs='+' list into the concrete list of axis values."""
    if len(vals) == 1:
        seq = [vals[0]]
    elif len(vals) == 3:
        start, step, stop = vals
        if step == 0:
            raise ValueError(f"--{name}: STEP must be non-zero for a 'START STEP STOP' sweep")
        eps = abs(step) * 1e-9
        n = int(math.floor((stop - start) / step + eps)) + 1
        if n < 1:
            raise ValueError(
                f"--{name}: STEP sign is inconsistent with START={start} STOP={stop}")
        seq = [start + i * step for i in range(n)]
    else:
        raise ValueError(
            f"--{name}: expected 1 value (fixed) or 3 values (START STEP STOP, "
            f"STOP inclusive); got {len(vals)}")
    if is_int:
        seq = [int(round(v)) for v in seq]
    return seq


def find_ltgs(explicit):
    """Locate the ltgs binary: explicit path, then ./build/ltgs next to the
    repo root, then $PATH."""
    if explicit:
        return explicit
    here = os.path.dirname(os.path.abspath(__file__))
    cand = os.path.join(here, os.pardir, "build", "ltgs")
    if os.path.exists(cand):
        return os.path.abspath(cand)
    from shutil import which
    found = which("ltgs")
    if found:
        return found
    raise SystemExit("could not find the ltgs binary; pass --ltgs-bin PATH")


def run_one(ltgs_bin, point):
    """Run ltgs for one parameter dict; return (E_per_spin, k_vec, k_idx, wts)."""
    with tempfile.TemporaryDirectory() as tmp:
        cmd = [ltgs_bin, str(int(point["L"]))]
        for name, _ in SCALAR_PARAMS:
            if name == "L" or name not in point:
                continue
            cmd += [f"--{name}", repr(point[name])]
        cmd += ["-o", tmp]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            raise RuntimeError(
                f"ltgs failed ({' '.join(cmd)}):\n{res.stderr}{res.stdout}")
        files = glob.glob(os.path.join(tmp, "*.h5"))
        if len(files) != 1:
            raise RuntimeError(f"expected exactly one output .h5, got {files}")
        with h5py.File(files[0], "r") as h:
            E_per_spin = float(h["E_min"][()])          # ltgs writes E/spin here
            k_vec = h["k_star_vec"][:].astype(float)
            k_idx = h["k_star_idx"][:].astype(np.int64)
            ev = h["eigenvector"][:]                     # (3*n_sl, 2) = (re, im)
            n_sl = h["sl_positions"].shape[0]
        comp = ev[:, 0] ** 2 + ev[:, 1] ** 2            # |component|^2
        wts = comp.reshape(n_sl, 3).sum(axis=1)         # per-sublattice weight
    return E_per_spin, k_vec, k_idx, wts


def main():
    ap = argparse.ArgumentParser(
        description="Sweep driver for the ltgs Luttinger-Tisza solver.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("L", nargs="+", type=float,
                    help="Linear supercell dimension (fixed, or START STEP STOP)")
    for name, _ in SCALAR_PARAMS:
        if name == "L":
            continue
        ap.add_argument(f"--{name}", nargs="+", type=float, default=None,
                        help="fixed value, or 'START STEP STOP' (STOP inclusive)")
    ap.add_argument("-o", "--output", required=True, help="output HDF5 file")
    ap.add_argument("--ltgs-bin", default=None, help="path to the ltgs binary")
    args = ap.parse_args()

    ltgs_bin = find_ltgs(args.ltgs_bin)

    # Expand each provided parameter into its concrete axis values.
    provided = {}   # name -> list of values, in SCALAR_PARAMS order
    for name, is_int in SCALAR_PARAMS:
        raw = args.L if name == "L" else getattr(args, name)
        if raw is None:
            continue
        provided[name] = expand(name, raw, is_int)

    swept = [n for n, _ in SCALAR_PARAMS if n in provided and len(provided[n]) > 1]
    fixed = [n for n, _ in SCALAR_PARAMS if n in provided and len(provided[n]) == 1]
    grid_shape = tuple(len(provided[n]) for n in swept)
    axis_values = [provided[n] for n in swept]
    fixed_point = {n: provided[n][0] for n in fixed}

    ntot = int(np.prod(grid_shape)) if grid_shape else 1
    print(f"ltgs: {ltgs_bin}", file=sys.stderr)
    print(f"swept axes {swept} -> grid {grid_shape or '(scalar)'}  "
          f"({ntot} run{'s' if ntot != 1 else ''}); fixed {fixed_point}",
          file=sys.stderr)

    # Iterate the Cartesian product in C order so a flat index reshapes cleanly.
    E = np.empty(ntot)
    lam = np.empty(ntot)
    kvec = np.empty((ntot, 3))
    kidx = np.empty((ntot, 3), dtype=np.int64)
    wts = None
    for flat, combo in enumerate(np.ndindex(*grid_shape) if grid_shape else [()]):
        point = dict(fixed_point)
        for ax, i in zip(swept, combo):
            point[ax] = provided[ax][i]
        E_i, kv_i, ki_i, w_i = run_one(ltgs_bin, point)
        if wts is None:
            n_sl = len(w_i)
            wts = np.empty((ntot, n_sl))
        E[flat] = E_i
        lam[flat] = 2.0 * E_i
        kvec[flat] = kv_i
        kidx[flat] = ki_i
        wts[flat] = w_i
        print(f"  [{flat + 1}/{ntot}] {point}  E/spin={E_i:.6f}", file=sys.stderr)

    # Reshape flat results onto the sweep grid.
    E = E.reshape(grid_shape)
    lam = lam.reshape(grid_shape)
    kvec = kvec.reshape(grid_shape + (3,))
    kidx = kidx.reshape(grid_shape + (3,))
    wts = wts.reshape(grid_shape + (n_sl,))

    with h5py.File(args.output, "w") as h:
        h.attrs["swept_params"] = ",".join(swept)
        h.attrs["grid_shape"] = np.array(grid_shape, dtype=np.int64)
        h.attrs["n_sl"] = n_sl
        h.attrs["ltgs_bin"] = ltgs_bin
        for name in fixed:
            h.attrs[f"fixed_{name}"] = provided[name][0]
        for name, vals in zip(swept, axis_values):
            h.create_dataset(f"axes/{name}", data=np.array(vals))
        h.create_dataset("E_per_spin", data=E)
        h.create_dataset("lambda_min", data=lam)
        h.create_dataset("k_star_vec", data=kvec)
        h.create_dataset("k_star_idx", data=kidx)
        h.create_dataset("sublattice_weights", data=wts)

    print(f"saved {ntot} result(s) to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
