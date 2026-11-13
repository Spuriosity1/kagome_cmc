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

Usage:
    plot_dw_order_params.py [FILE.out.h5 ...]
If no file is given, globs for *.out.h5 in this script's directory.
"""
import glob
import os
import sys

import h5py
import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))


def load(path):
    with h5py.File(path, "r") as h:
        g = h["domain_wall"]
        d = dict(
            x=g["x"][:].astype(float),
            Mz=g["M"][:, 2],
            chi=g["chi"][:],
            pin=int(g["pin_thickness"][()]),
            nx=g["x"].shape[0],
            D_sia=float(g["D_sia"][()]),
            Ddm=float(g["D"][()]),
            ddm=float(g["d"][()]),
            T=float(g["T_sample"][()]),
        )
        # Mask structurally-empty chirality slabs (odd x carry no triangle).
        cnt = g["chi_count"][:] if "chi_count" in g else np.ones(d["nx"])
        d["chi_mask"] = cnt > 0
    return d


def main(paths):
    if not paths:
        paths = sorted(glob.glob(os.path.join(HERE, "*.out.h5")))
    if not paths:
        sys.exit("no .out.h5 files found")

    fig, (ax_m, ax_c) = plt.subplots(2, 1, sharex=True, figsize=(7.5, 6.0))

    for p in paths:
        r = load(p)
        label = (rf"$D_\mathrm{{sia}}={r['D_sia']:g}$, "
                 rf"$D={r['Ddm']:g}$, $d={r['ddm']:g}$")
        ax_m.plot(r["x"], r["Mz"], "-", lw=1.6, label=label)
        m = r["chi_mask"]
        ax_c.plot(r["x"][m], r["chi"][m], ".-", lw=1.4, ms=5, label=label)

    # Shade the pinned Dirichlet end slabs and mark the nominal wall centre.
    r0 = load(paths[0])
    pin, nx = r0["pin"], r0["nx"]
    for ax in (ax_m, ax_c):
        ax.axvspan(-0.5, pin - 0.5, color="0.85", zorder=0)
        ax.axvspan(nx - pin - 0.5, nx - 0.5, color="0.85", zorder=0)
        ax.axvline(nx / 2 - 0.5, color="0.6", ls=":", lw=0.9)
        ax.axhline(0.0, color="0.6", lw=0.6)
        ax.margins(x=0)

    ax_m.set_ylabel(r"$\langle S_z \rangle$  (FM order parameter)")
    ax_c.set_ylabel(r"$\chi = \langle S_i\cdot(S_j\times S_k)\rangle$")
    ax_c.set_xlabel("x-slab index (wall normal)")
    ax_m.set_title("Domain-wall order parameters")
    if len(paths) > 1:
        ax_m.legend(fontsize=8, ncol=2)

    out = os.path.join(HERE, "dw_order_params.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("saved", out)


if __name__ == "__main__":
    main(sys.argv[1:])
