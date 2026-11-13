#!/usr/bin/env python3
"""Real-space domain-wall order parameters across a single-ion-anisotropy scan.

Plots the bg_ D-scan (Co3Sn2S2 parameters, Lx=20), one trace per run, against
the x-slab index (real space along the wall normal). Two stacked panels show the
two order parameters:

  top     <S_z>(x)   the FM magnetisation along the easy axis (reverses across
                     the wall), and
  bottom  chi(x)     the scalar spin chirality <S_i.(S_j x S_k)>, switched on by
                     the nearest-neighbour DMI and localised on the wall.

Trace colour encodes the single-ion anisotropy D_sia (log-scaled colormap).

Backward compatible with legacy files from the old `DC` binary: those store the
anisotropy under the key `D` (new files use `D_sia`) and carry no `chi` dataset,
so the chirality trace is simply skipped for runs that lack it.

Usage:
    plot_dw_profile.py [FILE.out.h5 ...]
With no arguments, globs for bg_*Lx=20*.out.h5 in this script's directory.
"""
import glob
import os
import sys

import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors

HERE = os.path.dirname(os.path.abspath(__file__))
PATTERN = os.path.join(HERE, "bg_*Lx=20*.out.h5")

# --- load ------------------------------------------------------------------
paths = sys.argv[1:] or sorted(glob.glob(PATTERN))
if not paths:
    sys.exit(f"no files matched {PATTERN!r} (and none given on the command line)")

runs = []
for f in paths:
    with h5py.File(f, "r") as h:
        g = h["domain_wall"]
        # Anisotropy key: new files store D_sia; legacy DC files store D.
        D_sia = float(g["D_sia"][()]) if "D_sia" in g else float(g["D"][()])
        rec = dict(
            D=D_sia,
            x=g["x"][:].astype(float),
            M=g["M"][:],                       # (nx, 3)
            pin=int(g["pin_thickness"][()]),
            nx=g["x"].shape[0],
            chi=None,
            chi_mask=None,
        )
        if "chi" in g:
            rec["chi"] = g["chi"][:]
            cnt = g["chi_count"][:] if "chi_count" in g else np.ones(rec["nx"])
            rec["chi_mask"] = cnt > 0    # structurally-empty slabs carry no triangle
        runs.append(rec)
runs.sort(key=lambda r: r["D"])
Ds = np.array([r["D"] for r in runs])
n_chi = sum(r["chi"] is not None for r in runs)
print(f"loaded {len(runs)} runs, D_sia = {Ds}  ({n_chi} with chirality)")

# --- colour map by D_sia (log scale) ---------------------------------------
norm = colors.LogNorm(vmin=Ds.min(), vmax=Ds.max())
cmap = cm.viridis
sm = cm.ScalarMappable(norm=norm, cmap=cmap)

# --- plot ------------------------------------------------------------------
fig, (ax_m, ax_c) = plt.subplots(2, 1, sharex=True, figsize=(7.5, 6.5))

pin = runs[0]["pin"]
nx = runs[0]["nx"]

# top panel: FM magnetisation <S_z>(x)
for r in runs:
    ax_m.plot(r["x"], r["M"][:, 2], "-", lw=1.6, color=cmap(norm(r["D"])))

# bottom panel: scalar spin chirality chi(x), masking empty slabs
for r in runs:
    if r["chi"] is None:
        continue
    m = r["chi_mask"]
    ax_c.plot(r["x"][m], r["chi"][m], "-", lw=1.4, color=cmap(norm(r["D"])))
if n_chi == 0:
    ax_c.text(0.5, 0.5, "no chi dataset in these files\n"
              "(re-run the domain_wall binary)",
              ha="center", va="center", transform=ax_c.transAxes, color="0.4")

for ax in (ax_m, ax_c):
    # shade pinned end slabs (Dirichlet BC), mark nominal wall centre + zero line
    ax.axvspan(-0.5, pin - 0.5, color="0.85", zorder=0)
    ax.axvspan(nx - pin - 0.5, nx - 0.5, color="0.85", zorder=0)
    ax.axvline(nx / 2 - 0.5, color="0.6", ls=":", lw=0.9)
    ax.axhline(0.0, color="0.6", lw=0.6)
    ax.margins(x=0)

ax_m.set_ylabel(r"$\langle S_z \rangle$  (FM order parameter)")
ax_c.set_ylabel(r"$\chi = \langle S_i\cdot(S_j\times S_k)\rangle$")
ax_c.set_xlabel("x-slab index (wall normal)")
ax_m.set_title("Domain-wall order parameters  (Co$_3$Sn$_2$S$_2$ params, $L_x=20$)")

cbar = fig.colorbar(sm, ax=(ax_m, ax_c), location="right", fraction=0.046, pad=0.02)
cbar.set_label("single-ion anisotropy $D_\\mathrm{sia}$")
cbar.set_ticks(Ds)
cbar.set_ticklabels([f"{d:g}" for d in Ds])

out = os.path.join(HERE, "dw_profile_bg.png")
fig.savefig(out, dpi=150, bbox_inches="tight")
print("saved", out)
