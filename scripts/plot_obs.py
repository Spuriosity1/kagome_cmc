#!/usr/bin/env python
"""Process a list of *.out.h5 run files (pass via shell glob) and plot a 3-panel
figure: heat capacity, FM order parameter (m_z), AFM order parameter (scalar
chirality).

Filename parameters (L, J1, Jc1, D, d, B, seed, Tc, sw, ...) are parsed from each
name. With -s PARAM [PARAM ...] the files are binned into one trace per distinct
combination of those parameters (e.g. -s Tc L), each labelled by its series name.
Files within a bin are blindly averaged -- no check that they otherwise match.

Put `--` before the file glob when using -s, so the parameter list and the files
don't run together:

usage: python plot_run2.py <file.out.h5> [...]
       python plot_run2.py -s Tc L -- out/L8/*.out.h5
"""
import argparse
import os
import re
from collections import defaultdict

import h5py
import numpy as np
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("-s", "--series", nargs="+", default=[],
                help="filename parameters to bin traces by, e.g. -s Tc L")
ap.add_argument("--save")
ap.add_argument("files", nargs="+", help="*.out.h5 run files (shell glob)")
args = ap.parse_args()


def parse_params(fname):
    """key=value pairs from a run filename (values are strings; B=0,0,0 kept whole).
    Keys/values are underscore-delimited, so neither may contain an underscore."""
    return dict(re.findall(r"([^_=]+)=([^_]+)", os.path.basename(fname)))


def _isnum(x):
    try:
        float(x)
        return True
    except ValueError:
        return False


def fmt_param(k, vals):
    """One value -> 'k=v'; several -> 'k in [lo, hi]' over the full range averaged."""
    vals = sorted(vals, key=lambda x: float(x) if _isnum(x) else x)
    if len(vals) == 1:
        return f"{k}={vals[0]}"
    return rf"{k}$\in$[{vals[0]}, {vals[-1]}]"


def load(fname):
    """Return (T, C, net_m_z, |chi|) for one run file, or None if it can't be read
    (e.g. a simulation still writing the file)."""
    L = int(parse_params(fname)["L"])
    N = 3 * L**3                                         # kagome: 3 sublattices per cell
    try:
        h = h5py.File(fname, "r")
    except OSError as e:
        print(f"skipping {os.path.basename(fname)}: {e}")
        return None
    with h:
        T = h["energy/T_list"][:]
        E = h["energy/E"][:]                             # sum of E over samples
        E2 = h["energy/E2"][:]                           # sum of E^2 over samples
        n = h["energy/n_samples"][:]
        C = N * (E2 / n - (E / n) ** 2) / T**2           # heat capacity
        mz = h["order/magnetization"][:, :, 2]           # (T, sublattice) z-component
        chi = h["order/scalar_chirality"][:]             # (T,)
    return T, C, np.abs(mz.mean(axis=1)), np.abs(chi)    # |net m_z|, |scalar chirality|


# bin files by the chosen series parameters (one bin if -s omitted)
bins = defaultdict(list)
for f in args.files:
    p = parse_params(f)
    key = tuple(f"{s}={p.get(s, '?')}" for s in args.series)
    bins[key].append(f)

fig, ax = plt.subplots(3, 1, sharex=True, figsize=(6, 8))
n_used = n_series = 0
param_vals = defaultdict(set)                             # filename params over used files

# colour each trace by the (first) series parameter's value; fall back to the
# default matplotlib colour cycle if any of those values won't parse as a float
keys = sorted(bins)
key_color = {}
if args.series:
    s0 = args.series[0]
    try:
        scalars = {key: float(dict(kv.split("=", 1) for kv in key)[s0]) for key in keys}
        print(scalars)
        min_s = min(scalars.values())
        max_s = max(scalars.values())
        if min_s * max_s < 0:
            cmap = plt.get_cmap("coolwarm")
            bound = max(max_s, -min_s)
            norm = plt.Normalize(-bound, bound)
        else:
            cmap = plt.get_cmap("viridis")
            norm = plt.Normalize(min_s, max_s)

        key_color = {key: cmap(norm(v)) for key, v in scalars.items()}
    except ValueError:
        key_color = {}

for key in keys:
    loaded = [(f, load(f)) for f in bins[key]]
    runs = [r for _, r in loaded if r is not None]
    if not runs:
        continue
    for f, r in loaded:
        if r is not None:
            for k, v in parse_params(f).items():
                param_vals[k].add(v)
    n_used += len(runs)
    n_series += 1
    T, C, m, chi = (np.mean(v, axis=0) for v in zip(*runs))
    label = ", ".join(key) if key else None
    color = key_color.get(key)                           # None -> default colour cycle
    ax[0].plot(T, C, ".-", label=label, color=color)
    ax[1].plot(T, m, ".-", label=label, color=color)
    ax[2].plot(T, chi, ".-", label=label, color=color)

ax[0].set_ylabel(r"Heat capacity $C$")
ax[1].set_ylabel(r"FM  $|m_z|$")
ax[2].set_ylabel(r"AFM  $|\chi_\mathrm{sc}|$")
ax[2].set_xlabel("Temperature $T$")
# describe the parameters folded into each trace (everything not used as a series key)
averaged = [fmt_param(k, param_vals[k])
            for k in sorted(param_vals) if k not in args.series]
title = f"{n_used} files, {n_series} series"
if averaged:
    title += "\n" + ", ".join(averaged)
ax[0].set_title(title)
for a in ax:
    a.grid(alpha=0.3)
if args.series:
    ax[0].legend(title=", ".join(args.series), fontsize="small")
fig.tight_layout()

if args.save is not None:
    fig.savefig(args.save)
    print("Saved figure to ",args.save)
else:
    plt.show()
