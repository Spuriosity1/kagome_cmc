# kagome_AC_CMC

Classical Monte Carlo for Heisenberg spins on the **kagome** lattice, using
[LatticeLab 2](https://github.com/Spuriosity1/latticelab) (`lattice_lib ≥ 2.3.0`).

Two protocols are provided:

- **`anneal`** — simulated annealing from `T_hot` down to `T_cold`, computing the
  static structure factor (SSF) S(q) via FFT at chosen sample temperatures.
- **`AC_sus`** — anneals to a target temperature, then drives the system with an
  oscillating field to measure the AC magnetic susceptibility χ′(ω) and χ″(ω).

Results are written to HDF5 for downstream analysis.

## Dependencies

| Library | Purpose |
|---------|---------|
| `lattice_lib ≥ 2.3.0` | Supercell geometry, periodic BCs, FFT infrastructure |
| `fftw3` | Underlying FFT engine |
| `hdf5` | Output file format (serial or parallel) |
| `nlohmann_json` | JSON (available for configuration) |
| `eigen3` | Linear algebra (used by the disabled `ltgs` target) |
| `meson ≥ 1.5.0` | Build system |
| Python ≥ 3.8 + `numpy`, `h5py`, `matplotlib` | Postprocessing / plotting |

`XoshiroCpp` (Xoshiro128++ PRNG) is bundled as a header in `include/`.

## Build

```bash
meson setup build
meson compile -C build
```

Outputs: `build/anneal` and `build/AC_sus` (both installed).

## Lattice

- **Kagome** — 3 sublattices per unit cell, at integer positions `(0,0,0)`,
  `(1,0,0)`, `(1,1,0)` (see `include/geometry/kagome.hpp`).
- The supercell is `[[2L,0,0],[L,2L,0],[0,0,L]]` in primitive-index units, giving
  **3·L³ spins** (3 sublattices × L³ kagome cells). The `L` dimension stacks
  layers along the out-of-plane (c) axis for the interlayer coupling `Jc1`.

### Couplings

All couplings are isotropic Heisenberg. Defaults are 0, so only the ones you pass
are activated.

| Flag | Neighbour |
|------|-----------|
| `--J1` | Nearest neighbour |
| `--J2` | Second neighbour |
| `--J3` | Third neighbour, across bonds |
| `--Jd` | Third neighbour, across hexagons |
| `--Jc1` | Out-of-plane (interlayer) |

**Hamiltonian:** H = Σ_n J_n Σ_{⟨ij⟩_n} S_i · S_j − B · Σ_i S_i (spins are unit
vectors on S²).

## Usage

### `anneal` — simulated annealing + SSF

Anneals from `T_hot` to `T_cold` in `n_steps` log-spaced steps. Energy statistics
are collected at every temperature; SSF statistics are collected only at the
`--T_sample` temperatures (which are inserted exactly into the cooling schedule).

```
anneal -o <output_dir> -s <seed> L <int>
       [--J1 <f>] [--J2 <f>] [--J3 <f>] [--Jd <f>] [--Jc1 <f>]
       --T_sample <float> [<float> ...]
       [--T_hot <float>] [--T_cold <float>] [--T_ref <float>]
       [--n_steps <int>] [--n_sweep <int>] [--n_burn_in <int>] [--n_sample <int>]
       [-B <bx> <by> <bz>] [--prefix <str>] [--save_state]
```

| Flag | Default | Description |
|------|---------|-------------|
| `-o / --output_dir` | required | Directory for output files (must exist) |
| `-s / --seed` | required | Seed index into the true-random cache (integer) |
| `L` | required | Supercell linear dimension; total spins = 3·L³ |
| `--J1 … --Jc1` | 0 | Exchange couplings (see table above) |
| `--T_sample` | required | One or more temperatures at which to collect SSF |
| `--T_hot` | 10 | Starting temperature |
| `--T_cold` | min `T_sample` | Final temperature |
| `--T_ref` | 1.0 | Reference temperature for Metropolis proposal width |
| `--n_steps` | 100 | Annealing temperature steps |
| `--n_sweep / -w` | 16 | Metropolis sweeps per temperature step |
| `--n_burn_in` | 16 | Burn-in sweeps at `T_hot` |
| `--n_sample` | 64 | Sampling rounds at each `T_sample` |
| `-B / --external_field` | (0,0,0) | External magnetic field vector |
| `--prefix` | `run` | Filename prefix |
| `--save_state` | — | Also write the final spin configuration to `.spins.h5` |

**Example:**
```bash
./build/anneal -o results/ -s 0 L 4 \
    --J1 1.0 --T_sample 0.05 --n_steps 200 --n_sample 128
```

**Output file:** `results/run_L=4_J1=1_B=0,0,0_seed=0_Tc=0.05_sw=16_.out.h5`
(fields are joined by `_`; only the couplings you pass appear in the name).

### `AC_sus` — AC susceptibility

Burns in at `T_hot`, anneals to `T_sample`, then applies an oscillating field
`B(t) = B_dc + sin(ω·n_sweep·i)·δB` over `n_steps` time steps (each step is
`n_sweep` MC sweeps), accumulating the magnetisation and decomposing it into
χ′(ω) and χ″(ω).

```
AC_sus -o <output_dir> -s <seed> L <int>
       [--J1 <f> … --Jc1 <f>]
       --T_sample <float> --omega <float>
       [--deltaB <dx> <dy> <dz>] [-B <bx> <by> <bz>]
       [--T_hot <float>] [--T_ref <float>] [--n_anneal_steps <int>]
       [--n_burn_in <int>] [--n_steps <int>] [--n_sweep <int>]
       [--prefix <str>] [--save_state]
```

| Flag | Default | Description |
|------|---------|-------------|
| `-o / --output_dir` | required | Directory for output files (must exist) |
| `-s / --seed` | required | Seed index into the true-random cache (integer) |
| `L`, `--J1 … --Jc1` | — / 0 | Lattice size and couplings (as for `anneal`) |
| `--T_sample` | required | Sampling temperature (final anneal target) |
| `--omega` | required | Drive frequency (rad / sweep) |
| `--deltaB` | (0,0,0) | AC field amplitude and direction |
| `-B / --external_field` | (0,0,0) | DC (bias) field |
| `--T_hot` | 10 | Temperature to begin annealing from |
| `--T_ref` | 1.0 | Reference temperature for Metropolis proposal width |
| `--n_anneal_steps` | 100 | Temperature grid points from `T_hot` to `T_sample` |
| `--n_burn_in` | 160 | Burn-in sweeps at `T_hot` |
| `--n_steps` | 1024 | AC time steps (each = `n_sweep` sweeps) |
| `--n_sweep / -w` | 16 | MC sweeps per AC time step |
| `--prefix` | `run` | Filename prefix |
| `--save_state` | — | Also write the final spin configuration to `.spins.h5` |

The susceptibility is extracted from the per-spin magnetisation M(i) as
```
χ′_α  =  (2/N) Σ_i  M_α(i) sin(phase_i) / |δB|
χ″_α  = −(2/N) Σ_i  M_α(i) cos(phase_i) / |δB|
M_dc_α =  (1/N) Σ_i  M_α(i)
```
with `phase_i = ω·n_sweep·i`. Positive χ″ corresponds to absorption.

**Example:**
```bash
./build/AC_sus -o results/ -s 0 L 4 \
    --J1 1.0 --T_sample 0.1 --omega 0.05 --deltaB 0.01 0 0 --n_steps 2048
```

## HDF5 output format

### `anneal` — `<params>.out.h5`

```
├── /geometry/
│   ├── recip_vectors[3, 3]     reciprocal lattice basis (rows: b0, b1, b2)
│   └── index_cell[3, 3]        supercell index matrix (int64)
├── /energy/
│   ├── T_list[n_T]
│   ├── E[n_T]                  sum of energy samples per temperature
│   ├── E2[n_T]                 sum of energy² samples
│   └── n_samples[n_T]
└── /ssf/
    ├── static_corr[n_corr, n_T, n_k, n_sl, n_sl, 2]   last dim: [re, im]
    ├── corr_lookup[n_corr]     e.g. ["xx", "yy", "zz"]
    ├── T_list[n_T]
    ├── n_samples[n_T]
    ├── sl_positions[n_sl, 3]   sublattice coordinates (int64)
    ├── @n_spins                total spin count (attribute)
    └── @k_dims[3]              k-space grid dimensions (attribute)
```

`static_corr` stores raw per-cell DFT correlators **without** sublattice phases;
these must be applied downstream (see [Postprocessing](#postprocessing)).

### `AC_sus` — `<params>.out.h5`

```
├── /ac_susceptibility/
│   ├── chi_prime[3]            in-phase susceptibility χ′ (per spin)
│   ├── chi_second[3]           out-of-phase susceptibility χ″ (per spin)
│   ├── M_dc[3]                 DC magnetisation (per spin)
│   ├── B[3], deltaB[3]         DC and AC field vectors
│   ├── deltaB_mag              |δB|
│   ├── T_sample, omega         scalars
│   └── n_steps, n_sweep        scalars (uint64)
├── /energy/                    annealing energy trace (as above)
└── /geometry/                  lattice metadata (as above)
```

## Postprocessing

The HDF5 outputs are self-describing and can be loaded with `h5py`. Two
conventions are needed when working with the SSF:

- **Sublattice phase correction.** `static_corr` stores raw per-cell DFT
  correlators without sublattice phases. To recover S(q), multiply by
  ```
  w[k, μ, ν] = exp(2πi Σ_j K_j (r_μ − r_ν)_j / k_dims_j)
  ```
  where `r_μ` are the sublattice positions (`/ssf/sl_positions`) and `k_dims`
  is the k-space grid (`@k_dims`), then sum over the sublattice indices μ, ν.
- **Γ centring.** `fftshift` the k-space grid so the Γ point appears at the
  centre.

Specific heat follows from the `/energy` group as
C_v = (⟨E²⟩ − ⟨E⟩²)/T².

## Code structure

```
include/
  MC.hpp                  MC_runner class: spins, bonds, Metropolis update
  cli_bits.hpp            CLI declarations, lattice/MC construction, T-profile
  ssf_manager.hpp         SSF accumulator: FFTs spin components, writes /ssf
  energy_manager.hpp      Energy accumulator, writes /energy
  abstract_manager.hpp    Base class for temperature-resolved accumulators
  geometry/kagome.hpp     Kagome geometry: sublattices and bond tables
  format_bits.hpp         Formatting utilities
  XoshiroCpp.hpp          Bundled Xoshiro128++ PRNG (header-only)
src/
  MC.cpp                  MC_runner implementation + HDF5 I/O
  anneal.cpp              Simulated annealing + SSF sampling
  AC_sus.cpp              AC susceptibility protocol
  luttinger_tisza.cpp     Luttinger–Tisza analysis (target currently disabled)
```

## Physics notes

- **Lattice:** kagome — 3 corner-sharing-triangle sublattices per cell, stacked
  along c for interlayer coupling.
- **Hamiltonian:** H = Σ_{⟨ij⟩_n} J_n S_i · S_j − B · Σ_i S_i (isotropic
  Heisenberg, in-plane up to third neighbours plus out-of-plane `Jc1`).
- **Metropolis proposal:** Gaussian kick in ℝ³ re-normalised to the unit sphere;
  width ∝ √(T/T_ref), with an over-relaxation step applied after each sweep.
- **SSF:** S^{αβ}(q) = (1/N) Σ_{μν} e^{iq·(r_μ−r_ν)} ⟨S^α_μ(q) S^{β*}_ν(q)⟩, with
  q on a uniform 3D grid set by `k_dims`.
```