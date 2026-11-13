#!/usr/bin/env python3
"""Spatial profiles of the two domain-wall order parameters.

Reads a `domain_wall` run (the binary renamed from DC) and plots, against the
x-slab index (the wall normal):

  OP1  <S_z>(x)   the FM magnetisation along the easy axis (the order parameter
                  that reverses across the wall), and
  OP2  chi(x)     the scalar spin chirality S_i.(S_j x S_k) averaged per slab,
                  which is switched on by the nearest-neighbour DMI and is
                  localised on the wall.

Chirality is binned by each triangle's cell-origin x (sublattice 0), which only
sits on even x; `chi_count==0` slabs carry no triangle and are masked out rather
than drawn as chi==0.

Each run carries its physical parameters (D_sia, D, d, T) in the file. Across a
set of files the parameters that are shared label the whole figure (title); the
one that varies is colour-coded (viridis, or coolwarm when it spans zero), the
same scheme as ../scripts/plot_obs.py. Use --color-by to force which parameter
drives the colour when several vary.

Usage:
    plot_dw_order_params.py [FILE.out.h5 ...]
    plot_dw_order_params.py --color-by T -- run_*.out.h5
If no file is given, globs for *.out.h5 in this script's directory.
"""
import argparse
import glob
import os
import sys
from collections import defaultdict

import h5py
import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))

# physical parameters stored per run, in display order: param key -> latex label.
# The HDF5 dataset name equals the key except where DSET overrides it.
PARAMS = {
    "J1": r"J_1",
    "J2": r"J_2",
    "J3": r"J_3",
    "Jd": r"J_d",
    "Jc1": r"J_\mathrm{c1}",
    "D_sia": r"D_\mathrm{sia}",
    "D": r"D",
    "d": r"d",
    "T": r"T",
}
DSET = {"T": "T_sample"}                                 # param key -> dataset name
COUPLINGS = {"J1", "J2", "J3", "Jd", "Jc1"}              # shown in title only if "live"


def load(path):
    with h5py.File(path, "r") as h:
        g = h["domain_wall"]
        # Missing datasets read as 0.0 so files written before a parameter was
        # added (e.g. older runs without the exchange couplings) still load.
        def param(k):
            name = DSET.get(k, k)
            return float(g[name][()]) if name in g else 0.0
        d = dict(
            x=g["x"][:].astype(float),
            Mz=g["M"][:, 2],
            chi=g["chi"][:],
            pin=int(g["pin_thickness"][()]),
            nx=g["x"].shape[0],
            params={k: param(k) for k in PARAMS},
        )
        # Mask structurally-empty chirality slabs (odd x carry no triangle).
        cnt = g["chi_count"][:] if "chi_count" in g else np.ones(d["nx"])
        d["chi_mask"] = cnt > 0
    return d


def fmt(k, v):
    """'D_sia=0.1'-style latex fragment for parameter k with value v."""
    return rf"${PARAMS[k]}={v:g}$"


def average_seeds(members):
    """Collapse runs sharing all physical parameters (typically different seeds)
    into one by averaging the per-slab profiles. Arrays are assumed aligned --
    nonsense if the geometry (x grid) actually differs between members."""
    r = dict(members[0])                                 # keep params, pin, nx, x
    r["Mz"] = np.mean([m["Mz"] for m in members], axis=0)
    r["chi"] = np.mean([m["chi"] for m in members], axis=0)
    r["chi_mask"] = np.logical_and.reduce([m["chi_mask"] for m in members])
    r["n_seeds"] = len(members)
    return r


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("-c", "--color-by", choices=list(PARAMS),
                    help="parameter to colour traces by (default: the one that varies)")
    ap.add_argument("--save")
    ap.add_argument("files", nargs="*", help="*.out.h5 run files (shell glob)")
    args = ap.parse_args(argv)

    paths = args.files or sorted(glob.glob(os.path.join(HERE, "*.out.h5")))
    if not paths:
        sys.exit("no .out.h5 files found")

    runs = [load(p) for p in paths]

    # average replicate runs that share every physical parameter (seed is not
    # among them, so differing-seed files bin together) into one trace
    groups = defaultdict(list)
    for r in runs:
        groups[tuple(r["params"][k] for k in PARAMS)].append(r)
    runs = [average_seeds(members) for members in groups.values()]
    n_files = sum(r["n_seeds"] for r in runs)
    print(f"{n_files} files -> {len(runs)} averaged traces")

    # split parameters into those shared across all files (-> title) and those
    # that vary (-> distinguish/colour the traces)
    vals = defaultdict(set)
    for r in runs:
        for k, v in r["params"].items():
            vals[k].add(f"{v:g}")
    varying = [k for k in PARAMS if len(vals[k]) > 1]
    # title gets the shared params, but drop exchange couplings that are zero
    # everywhere -- those are "dead" bonds, not worth the clutter
    common = [k for k in PARAMS if k not in varying
              and not (k in COUPLINGS and runs[0]["params"][k] == 0)]

    # colour by the chosen (or first varying) parameter, viridis for a one-sided
    # range and coolwarm when it straddles zero; constant param -> default cycle
    ckey = args.color_by or (varying[0] if varying else None)
    color_of = {}
    if ckey is not None:
        s = [r["params"][ckey] for r in runs]
        mn, mx = min(s), max(s)
        if mn * mx < 0:
            cmap, bound = plt.get_cmap("coolwarm"), max(mx, -mn)
            norm = plt.Normalize(-bound, bound)
        elif mx > mn:
            cmap, norm = plt.get_cmap("jet"), plt.Normalize(mn, mx)
        else:
            cmap = norm = None
        if cmap is not None:
            color_of = {id(r): cmap(norm(r["params"][ckey])) for r in runs}

    fig, (ax_m, ax_c) = plt.subplots(2, 1, sharex=True, figsize=(7.5, 6.0))

    # draw ordered by the colour parameter so the legend reads monotonically
    order = sorted(runs, key=lambda r: r["params"][ckey]) if ckey else runs
    for r in order:
        label = ", ".join(fmt(k, r["params"][k]) for k in varying) or None
        color = color_of.get(id(r))
        ax_m.plot(r["x"], r["Mz"], "-", lw=1.6, color=color, label=label)
        m = r["chi_mask"]
        ax_c.plot(r["x"][m], r["chi"][m], ".-", lw=1.4, ms=5, color=color, label=label)

    # Shade the pinned Dirichlet end slabs and mark the nominal wall centre.
    pin, nx = runs[0]["pin"], runs[0]["nx"]
    for ax in (ax_m, ax_c):
        ax.axvspan(-0.5, pin - 0.5, color="0.85", zorder=0)
        ax.axvspan(nx - pin - 0.5, nx - 0.5, color="0.85", zorder=0)
        ax.axvline(nx / 2 - 0.5, color="0.6", ls=":", lw=0.9)
        ax.axhline(0.0, color="0.6", lw=0.6)
        ax.margins(x=0)

    ax_m.set_ylabel(r"$\langle S_z \rangle$  (FM order parameter)")
    ax_c.set_ylabel(r"$\chi = \langle S_i\cdot(S_j\times S_k)\rangle$")
    ax_c.set_xlabel("x-slab index (wall normal)")

    # title: the parameters shared by every file on the figure
    title = "Domain-wall order parameters"
    shared = ", ".join(fmt(k, runs[0]["params"][k]) for k in common)
    if shared:
        title += "\n" + shared
    ax_m.set_title(title)
    if varying:
        ax_m.legend(title=", ".join(f"${PARAMS[k]}$" for k in varying),
                    fontsize=8, ncol=2)

    out = args.save or os.path.join(HERE, "dw_order_params.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("saved", out)
    if args.save is None:
        plt.show()


if __name__ == "__main__":
    main(sys.argv[1:])
