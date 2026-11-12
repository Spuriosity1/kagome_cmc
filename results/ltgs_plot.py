#!/usr/bin/env python3
"""Plot the output of `ltgs_sweep.py`.

Reads the self-describing HDF5 file written by the sweep driver and renders the
Luttinger-Tisza results against the swept parameter(s):

  * 0 swept axes  -> prints a one-line summary (nothing to plot).
  * 1 swept axis  -> stacked line panels: E/spin, k* components, sublattice
                     weights, all vs the swept parameter.
  * 2 swept axes  -> heatmaps over the two axes: E/spin, |k*|, dominant
                     sublattice weight.
  * >2 swept axes -> unsupported (slice the file first).

Usage:
    ./ltgs_plot.py sweep.h5 [-o figure.png]
"""
import argparse
import os

import h5py
import numpy as np
import matplotlib.pyplot as plt


def _s(x):
    """Decode an HDF5 attribute that may come back as bytes."""
    return x.decode() if isinstance(x, bytes) else str(x)


def load(path):
    with h5py.File(path, "r") as h:
        swept = [s for s in _s(h.attrs["swept_params"]).split(",") if s]
        data = dict(
            swept=swept,
            n_sl=int(h.attrs["n_sl"]),
            axes={n: h[f"axes/{n}"][:] for n in swept},
            E=h["E_per_spin"][()],
            lam=h["lambda_min"][()],
            kvec=h["k_star_vec"][()],
            kidx=h["k_star_idx"][()],
            wts=h["sublattice_weights"][()],
            fixed={k[len("fixed_"):]: h.attrs[k]
                   for k in h.attrs if k.startswith("fixed_")},
        )
    return data


def fixed_title(d):
    bits = [f"{k}={v:g}" for k, v in d["fixed"].items()]
    return "ltgs sweep   " + ", ".join(bits) if bits else "ltgs sweep"


def plot_1d(d, out):
    name = d["swept"][0]
    x = d["axes"][name]
    kvec = np.atleast_2d(d["kvec"])        # (nx, 3)
    wts = np.atleast_2d(d["wts"])          # (nx, n_sl)

    fig, (ax0, ax1, ax2) = plt.subplots(3, 1, sharex=True, figsize=(7.5, 8.5))

    ax0.plot(x, d["E"], "-o", lw=1.6, ms=4, color="C0")
    ax0.set_ylabel(r"$E/\mathrm{spin}$")
    ax0.axhline(0.0, color="0.6", lw=0.6)

    for c, lbl in enumerate(("$k^*_x$", "$k^*_y$", "$k^*_z$")):
        ax1.plot(x, kvec[:, c], "-o", lw=1.4, ms=3, label=lbl)
    ax1.set_ylabel(r"$k^*$  (rad/coord-unit)")
    ax1.legend(ncol=3, fontsize=9)

    for s in range(wts.shape[1]):
        ax2.plot(x, wts[:, s], "-o", lw=1.4, ms=3, label=f"SL {s}")
    ax2.set_ylabel(r"$|s_\alpha|^2$ at $k^*$")
    ax2.legend(ncol=wts.shape[1], fontsize=9)
    ax2.set_xlabel(name)

    for ax in (ax0, ax1, ax2):
        ax.margins(x=0.02)
    ax0.set_title(fixed_title(d))
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("saved", out)


def _heat(ax, xname, yname, xv, yv, C, title, cmap="viridis"):
    pc = ax.pcolormesh(xv, yv, C.T, shading="auto", cmap=cmap)
    ax.figure.colorbar(pc, ax=ax, fraction=0.046, pad=0.03)
    ax.set_xlabel(xname)
    ax.set_ylabel(yname)
    ax.set_title(title)


def plot_2d(d, out):
    xn, yn = d["swept"]
    xv, yv = d["axes"][xn], d["axes"][yn]
    kmag = np.linalg.norm(d["kvec"], axis=-1)     # (nx, ny)
    wmax = d["wts"].max(axis=-1)                   # dominant sublattice weight

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    _heat(axes[0], xn, yn, xv, yv, d["E"], r"$E/\mathrm{spin}$")
    _heat(axes[1], xn, yn, xv, yv, kmag, r"$|k^*|$", cmap="magma")
    _heat(axes[2], xn, yn, xv, yv, wmax, r"$\max_\alpha |s_\alpha|^2$",
          cmap="cividis")
    fig.suptitle(fixed_title(d))
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("saved", out)


def summarize(d):
    print("no swept axes; single ltgs result:")
    print(f"  fixed      : {d['fixed']}")
    print(f"  E/spin     : {float(d['E']):.6f}")
    print(f"  lambda_min : {float(d['lam']):.6f}")
    print(f"  k*_vec     : {np.asarray(d['kvec'])}")
    print(f"  k*_idx     : {np.asarray(d['kidx'])}")
    print(f"  |s_a|^2    : {np.asarray(d['wts'])}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sweep", help="HDF5 file written by ltgs_sweep.py")
    ap.add_argument("-o", "--output", default=None,
                    help="output image (default: <sweep>.png)")
    args = ap.parse_args()

    d = load(args.sweep)
    out = args.output or os.path.splitext(args.sweep)[0] + ".png"

    ndim = len(d["swept"])
    if ndim == 0:
        summarize(d)
    elif ndim == 1:
        plot_1d(d, out)
    elif ndim == 2:
        plot_2d(d, out)
    else:
        raise SystemExit(
            f"{ndim}-D sweep over {d['swept']} is not plottable directly; "
            "slice the HDF5 file down to 1 or 2 axes first.")


if __name__ == "__main__":
    main()
