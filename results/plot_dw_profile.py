#!/usr/bin/env python3
"""Three-panel domain-wall magnetisation profile <S>(x), one trace per D value.

Plots the bg_ D-scan (Co3Sn2S2 parameters, Lx=20). Panels = Mx, My, Mz;
trace colour encodes the single-ion anisotropy D (log-scaled colormap).
"""
import glob
import os

import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors

HERE = os.path.dirname(os.path.abspath(__file__))
PATTERN = os.path.join(HERE, "bg_*Lx=20*.out.h5")

# --- load ------------------------------------------------------------------
runs = []
for f in glob.glob(PATTERN):
    with h5py.File(f, "r") as h:
        g = h["domain_wall"]
        runs.append(dict(
            D=float(g["D"][()]),
            x=g["x"][:].astype(float),
            M=g["M"][:],                       # (nx, 3)
            pin=int(g["pin_thickness"][()]),
            nx=g["x"].shape[0],
        ))
runs.sort(key=lambda r: r["D"])
Ds = np.array([r["D"] for r in runs])
print(f"loaded {len(runs)} runs, D = {Ds}")

# --- colour map by D (log scale, spans 1..60) ------------------------------
norm = colors.LogNorm(vmin=Ds.min(), vmax=Ds.max())
cmap = cm.viridis
sm = cm.ScalarMappable(norm=norm, cmap=cmap)

# --- plot ------------------------------------------------------------------
comps = ["M_x", "M_y", r"M_z \;(\mathrm{order\ parameter})"]
fig, axes = plt.subplots(3, 1, sharex=True, figsize=(7.5, 8.5))

pin = runs[0]["pin"]
nx = runs[0]["nx"]
for k, ax in enumerate(axes):
    for r in runs:
        ax.plot(r["x"], r["M"][:, k], "-", lw=1.6,
                color=cmap(norm(r["D"])))
    # shade pinned end slabs (Dirichlet BC)
    ax.axvspan(-0.5, pin - 0.5, color="0.85", zorder=0)
    ax.axvspan(nx - pin - 0.5, nx - 0.5, color="0.85", zorder=0)
    ax.axvline(nx / 2 - 0.5, color="0.6", ls=":", lw=0.9)
    ax.axhline(0.0, color="0.6", lw=0.6)
    ax.set_ylabel(f"$\\langle {comps[k]} \\rangle$")
    ax.margins(x=0)

axes[-1].set_xlabel("x-slab index")
axes[0].set_title("Domain-wall magnetisation profile  (Co$_3$Sn$_2$S$_2$ params, $L_x=20$)")

cbar = fig.colorbar(sm, ax=axes, location="right", fraction=0.046, pad=0.02)
cbar.set_label("single-ion anisotropy $D$")
cbar.set_ticks(Ds)
cbar.set_ticklabels([f"{d:g}" for d in Ds])

out = os.path.join(HERE, "dw_profile_bg.png")
fig.savefig(out, dpi=150, bbox_inches="tight")
print("saved", out)
