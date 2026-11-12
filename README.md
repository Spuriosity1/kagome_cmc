# kagome_AC_CMC

Classical Monte Carlo for Heisenberg spins on the **kagome** lattice, using
[LatticeLab 2](https://github.com/Spuriosity1/latticelab) (`lattice_lib ≥ 2.3.0`).

The following executables are provided:

- **`anneal`** — simulated annealing from `T_hot` down to `T_cold`, computing the
  static structure factor (SSF) S(q) via FFT at chosen sample temperatures.
- **`AC_sus`** — anneals to a target temperature, then drives the system with an
  oscillating field to measure the AC magnetic susceptibility χ′(ω) and χ″(ω).
- **`DC`** — static-field protocol for **domain walls**: builds an orthogonal
  cuboid elongated along the wall normal, pins the two end slabs to opposite
  easy-axis directions, and measures the equilibrium magnetisation profile. Adds
  a by-hand single-ion easy-axis anisotropy so the discrete wall can be driven
  into the Bulaevskii–Ginzburg (narrow-wall) regime. Parameterised for Co₃Sn₂S₂.
- **`ltgs`** — zero-temperature **Luttinger–Tisza** solver: sweeps the Brillouin
  zone, diagonalising the Fourier-transformed coupling matrix to find the ordering
  wavevector `k*` and the LT lower-bound energy per spin. No Monte Carlo is run.

Results are written to HDF5 for downstream analysis.

## Dependencies

| Library | Purpose |
|---------|---------|
| `lattice_lib ≥ 2.3.0` | Supercell geometry, periodic BCs, FFT infrastructure |
| `fftw3` | Underlying FFT engine |
| `hdf5` | Output file format (serial or parallel) |
| `nlohmann_json` | JSON (available for configuration) |
| `eigen3` | Linear algebra (BZ diagonalisation in `ltgs`) |
| `meson ≥ 1.5.0` | Build system |
| Python ≥ 3.8 + `numpy`, `h5py`, `matplotlib` | Postprocessing / plotting |

`XoshiroCpp` (Xoshiro128++ PRNG) is bundled as a header in `include/`.

## Build

```bash
meson setup build
meson compile -C build
```

Outputs: `build/anneal`, `build/AC_sus`, `build/DC`, and `build/ltgs` (all installed).

## Lattice

- **Kagome** — 3 sublattices per unit cell, at integer positions `(0,0,0)`,
  `(1,0,0)`, `(1,1,0)` (see `include/geometry/kagome.hpp`).
- The supercell is `[[2L,0,0],[L,2L,0],[0,0,L]]` in primitive-index units, giving
  **3·L³ spins** (3 sublattices × L³ kagome cells). The `L` dimension stacks
  layers along the out-of-plane (c) axis for the interlayer coupling `Jc1`.
- Stacking follows Co₃Sn₂S₂ (rhombohedral R‑3m): layers are ABC‑stacked, so each
  `Jc1` bond connects *different* sublattices (two neighbours per site in the layer
  above, two below). Only the interlayer bond topology is physical — the integer
  embedding leaves every layer at the same in-plane coordinates (a shear of the
  true ABC structure).

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

The burn-in and annealing stages use overrelaxation for fast equilibration.
The AC drive stage, by contrast, uses **pure Metropolis with a fixed proposal
width** (`--step_size`) and no overrelaxation, so that one sweep corresponds to
a fixed increment of (overdamped, Model-A) physical time — the MC-time ≈
real-time approximation appropriate for the kHz timescales of an AC experiment,
where spin precession is fully motionally averaged.

```
AC_sus -o <output_dir> -s <seed> L <int>
       [--J1 <f> … --Jc1 <f>]
       --T_sample <float> --omega <float>
       [--deltaB <dx> <dy> <dz>] [-B <bx> <by> <bz>]
       [--T_hot <float>] [--T_ref <float>] [--n_anneal_steps <int>]
       [--n_burn_in <int>] [--n_steps <int>] [--n_sweep <int>] [--step_size <f>]
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
| `--step_size` | 0.1 | Metropolis proposal width for the drive stage; sets the physical time per sweep (smaller = finer resolution, slower evolution) |
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

### `DC` — static field / domain wall

Builds an **orthogonal cuboid** supercell (independent `--Lx/--Ly/--Lz`, each
defaulting to the positional `L`) elongated along **x**, the wall normal. The
easy axis is **c (z)**. The two end slabs in x (thickness `--pin_thickness`
integer units) are **pinned** to +z and −z, imposing Dirichlet boundary
conditions that stabilise a single domain wall; y and z stay periodic. Because
the maximum bond range in x is 2 integer units, a pinned slab of thickness ≥ 2
makes the bulk behave as open-BC (wrap-around bonds across the x-seam only ever
connect frozen spins — no bond pruning needed).

A by-hand single-ion easy-axis anisotropy is added:
```
H = Σ_{<ij>_n} J_n S_i·S_j  −  B·Σ_i S_i  −  D Σ_i (S_i·ê)²     (ê = --easy_axis)
```
`D>0` is easy-axis. All stages use **pure Metropolis** (`sweep_metropolis`, no
overrelaxation), which is required for correct sampling once the on-site
anisotropy term is present. The protocol burns in at `T_hot`, anneals to
`T_sample` in the static field `B`, then averages the x-resolved profile ⟨S⟩(x)
over `n_sample` measurements.

```
DC -o <output_dir> -s <seed> L <int>
   [--Lx <int>] [--Ly <int>] [--Lz <int>]
   [--J1 <f> … --Jc1 <f>] [-B <bx> <by> <bz>]
   --T_sample <float>
   [--D <float>] [--easy_axis <ex> <ey> <ez>] [--pin_thickness <int>]
   [--T_hot <float>] [--T_ref <float>] [--n_anneal_steps <int>]
   [--n_burn_in <int>] [--n_sweep <int>] [--n_sample <int>]
   [--prefix <str>] [--save_state]
```

| Flag | Default | Description |
|------|---------|-------------|
| `-o / --output_dir` | required | Directory for output files (must exist) |
| `-s / --seed` | required | Seed index into the true-random cache (integer) |
| `L`, `--J1 … --Jc1` | — / 0 | Fallback size and couplings (as for `anneal`) |
| `--Lx / --Ly / --Lz` | `L` | Independent cuboid dimensions (cells); elongate `Lx` |
| `--T_sample` | required | Target temperature for the profile measurement |
| `--D` | 1.36 | Single-ion easy-axis anisotropy strength (meV) |
| `--easy_axis` | (0,0,1) | Easy-axis direction (need not be normalised) |
| `--pin_thickness` | 2 | Thickness (integer x units) of each pinned end slab (≥ 2) |
| `-B / --external_field` | (0,0,0) | Static (DC) field |
| `--T_hot` | 40 | Temperature to begin annealing from |
| `--T_ref` | 1.0 | Reference temperature for Metropolis proposal width |
| `--n_anneal_steps` | 200 | Temperature grid points from `T_hot` to `T_sample` |
| `--n_burn_in` | 200 | Burn-in sweeps at `T_hot` |
| `--n_sweep / -w` | 16 | Sweeps per anneal step / between profile samples |
| `--n_sample` | 200 | Profile measurements to average at `T_sample` |
| `--prefix` | `run` | Filename prefix |
| `--save_state` | — | Also write the full 3-D spin configuration to `.spins.h5` |

Units are meV (k_B = 1); Co₃Sn₂S₂ has T_C ≈ 175 K ≈ 15 meV.

**Co₃Sn₂S₂ parameters.** Exchange constants are taken from Zhang *et al.*,
[PRL **127**, 117201 (2021)](https://doi.org/10.1103/PhysRevLett.127.117201),
Table I (convention `H = Σ J_ij S_i·S_j`, **negative = FM**; values are `SJ` in
meV, used directly for unit spins). The documented default is the **MFWF
dominant-`J1`** set (Ref. [23] therein), which is a genuine collinear
ferromagnet in this classical model and therefore supports a clean FM domain
wall:

```bash
# MFWF dominant-J1 FM set — moderate D (near the BG crossover)
./build/DC -o results/ -s 0 4 --Lx 24 --Ly 3 --Lz 3 \
    --J1 -13.94 --J2 -0.94 --J3 -0.12 --Jd -0.12 --Jc1 -1.76 \
    --D 5 --T_hot 40 --T_sample 3 --save_state
```

> **Caveat on the "Fits to SW spectra" set.** The paper's headline set
> (`J2=+5.78, Jd=−29.94, Jc1=−11.56`, with `J1≈0`) does **not** stabilise the
> collinear ferromagnet in a classical model: the dominant `Jd` (across-hexagon)
> and `Jc1` are *intra-sublattice*, so the only inter-sublattice coupling left is
> the **AFM `J2`**, which anti-aligns the three kagome sublattices. The result is
> an even/odd-x staggered pattern rather than two FM domains — unsuitable for a
> clean FM domain-wall study. Use a dominant-`J1` set for domain walls.

**Bulaevskii–Ginzburg regime.** For an easy-axis ferromagnet the wall width
scales as `δ ~ a·√(J_stiff/D)`, so the narrow (Ising) wall is reached when
`D ≳ J_stiff` (the exchange stiffness along the wall normal). Scan `--D` (e.g.
`1, 5, 10, 20, 40` meV) and watch the ⟨S_z⟩(x) transition collapse from many
slabs to ~1–2 slabs.

**Limitation.** The kagome bond model has no second-interlayer `Jc2`; the
SW-fit set's small `Jc2=−1.70 meV` is omitted. Adding a `{{0,0,2}}` entry to
`kagome::vecs` plus a `--Jc2` flag is a one-line extension if required.

### `ltgs` — Luttinger–Tisza analysis

Builds the Fourier-transformed coupling matrix `J(k)` — an `n_sl × n_sl` complex
Hermitian matrix in sublattice space (the three spin components are degenerate
since all couplings are isotropic Heisenberg) — on the `L³` k-grid of a diagonal
supercell, and diagonalises it at every `k`. The global minimum eigenvalue `λ_min`
gives the Luttinger–Tisza lower bound on the energy, `E/spin = ½·λ_min` (the ½
accounts for each undirected bond appearing twice in the directed-bond sweep); the
`k` that attains it is the candidate ordering wavevector `k*`, and the associated
eigenvector gives the per-sublattice amplitude/phase pattern. No Monte Carlo is
run and no supercell of spins is allocated. The energy is directly comparable to
the per-spin energy from `anneal` (same bond tables, same convention).

The ordering-wavevector bound is only a true ground state if the LT eigenvector
can be realised with *every* spin a unit vector (the "strong" LT constraint);
otherwise `E/spin` is a lower bound. The couplings use the same flags and bond
tables as the MC protocols.

```
ltgs L <int>
     [--J1 <f>] [--J2 <f>] [--J3 <f>] [--Jd <f>] [--Jc1 <f>]
     [-o <output_dir>]
```

| Flag | Default | Description |
|------|---------|-------------|
| `L` | required | k-grid is `L³` (diagonal supercell in primitive-cell units) |
| `--J1 … --Jc1` | 0 | Exchange couplings (as for `anneal`) |
| `-o / --output_dir` | — | If given (must exist), also writes an HDF5 file; otherwise prints only |

`k*`, `λ_min`, and the per-sublattice eigenvector weights are printed to stdout.

**Example:**
```bash
./build/ltgs 16 --J1 1.0 --J2 0.3 --Jc1 -0.5 -o results/
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

### `DC` — `<params>.out.h5`

```
├── /domain_wall/
│   ├── x[nx]                   integer x-slab coordinate (nx = x-extent)
│   ├── M[nx, 3]                slab-averaged ⟨S⟩ profile (column 2 = M_z = OP)
│   ├── M_net[3]                net magnetisation (per spin)
│   ├── easy_axis[3], B[3]      easy-axis direction and static field
│   ├── D, T_sample             scalars (double)
│   └── Lx, Ly, Lz, pin_thickness   scalars (int)
├── /energy/                    annealing energy trace (as above)
└── /geometry/                  lattice metadata (as above)
```

The domain wall is read off the `M_z = M[:,2]` column as a function of `x`; its
width (and hence the Bulaevskii–Ginzburg regime) follows from a `tanh` fit.

### `ltgs` — `<params>lt.out.h5`

```
├── eigenvalue_map[Nk]         λ_min(k) over the flattened L³ k-grid
├── E_min                      scalar: LT energy per spin (= ½·min_k λ_min)
├── k_star_idx[3]              integer grid index of the minimising k (int64)
├── k_star_vec[3]              k* in reciprocal coordinates (rad / coord-unit)
├── eigenvector[n_sl, 2]       minimising eigenvector; last dim [re, im]
├── sl_positions[n_sl, 3]      sublattice coordinates (int64)
├── /geometry/
│   ├── index_cell[3, 3]       supercell index matrix (int64)
│   └── recip_vectors[3, 3]    reciprocal lattice basis (rows: b0, b1, b2)
└── @k_dims[3]                 k-space grid dimensions (attribute)
```

`eigenvalue_map` is flattened in the same row-major order as the SSF grids, so it
`reshape`s to `@k_dims` and `fftshift`es to centre Γ. Written only when `-o` is
given.

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
The raw correlators are stored on their own (as opposed to computing C on the fly)
to allow for Bessel corrections to be applied.

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
  DC.cpp                  Static-field domain-wall protocol (cuboid + anisotropy)
  luttinger_tisza.cpp     Luttinger–Tisza BZ sweep (ltgs); no MC, Eigen-based
```

## Physics notes

- **Lattice:** kagome — 3 corner-sharing-triangle sublattices per cell, stacked
  along c for interlayer coupling.
- **Hamiltonian:** H = Σ_{⟨ij⟩_n} J_n S_i · S_j − B · Σ_i S_i (isotropic
  Heisenberg, in-plane up to third neighbours plus out-of-plane `Jc1`). The `DC`
  protocol additionally adds a single-ion easy-axis term −D Σ_i (S_i·ê)².
- **Anisotropy & sampling:** with a single-ion term present, over-relaxation is
  no longer microcanonical, so `DC` uses pure-Metropolis equilibration
  (`sweep_metropolis`). Spins flagged `pinned` are frozen (Dirichlet BCs).
- **Metropolis proposal:** Gaussian kick in ℝ³ re-normalised to the unit sphere.
  Equilibration sweeps (`sweep_equilibrate`) use an adaptive width ∝ √(T/T_ref)
  followed by a microcanonical over-relaxation pass for fast decorrelation.
  Dynamics sweeps (`sweep_dynamics`, used by `AC_sus`'s drive stage) use a fixed
  proposal width and **no** over-relaxation, so each sweep is a fixed increment
  of overdamped physical time.
- **SSF:** S^{αβ}(q) = (1/N) Σ_{μν} e^{iq·(r_μ−r_ν)} ⟨S^α_μ(q) S^{β*}_ν(q)⟩, with
  q on a uniform 3D grid set by `k_dims`.
```
