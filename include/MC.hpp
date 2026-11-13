#pragma once
#include "H5Ipublic.h"
#include "H5Ppublic.h"
#include "H5Tpublic.h"
#include <nlohmann/json.hpp>
#include <filesystem>
#include <map>
#include <vector>
#include <string>
#include <cassert>

// LatticeLab 2
#include <lattice_lib/supercell.hpp>
#include <XoshiroCpp.hpp>
#include <random>


namespace CMC {

using json=nlohmann::json;


struct HeisenbergSpin;

/**
 * @brief Neighbours of a spin for one coupling shell, each with its own
 *        (already oriented) coupling matrix.
 *
 * The contribution of this shell to the spin's local field is
 *     h += Σ_i coupling[i] · neighbours[i]->S.
 * For an oriented bond source→target carrying matrix M (H_bond = S_srcᵀ M S_tgt)
 * the source stores M for the target, and the target stores Mᵀ for the source,
 * so each undirected bond is felt — with the correct orientation — by both ends.
 */
struct NeighbourSpins {
    std::vector<HeisenbergSpin*>             neighbours;
    std::vector<vector3::mat33<double>>      coupling;
};


/**
 * @brief Classical Heisenberg spin stored in the supercell.
 *
 * Satisfies the GeometricObject concept via the `ipos` member.
 */
struct HeisenbergSpin {
    ipos_t ipos;                           // required by GeometricObject concept
    int pyro_sl = 0;                       // kagome sublattice 0–2
    vector3::vec3<double> S;
    std::vector<NeighbourSpins> bond_sets;
    // If true this spin is frozen: Metropolis/overrelaxation skip it. Used to
    // impose Dirichlet boundary conditions that pin a domain wall (see domain_wall.cpp).
    bool pinned = false;
};

/**
 * @brief Supercell of Heisenberg spins on the kagome lattice.
 */
typedef Supercell<HeisenbergSpin> Lattice;


/**
 * @brief Specification of a single coupling term in the Hamiltonian.
 *
 * relative_vectors[sl] lists displacement vectors from a spin of sublattice sl
 * to its coupled neighbours (each undirected bond listed once; see kagome.hpp).
 *
 * The coupling matrix of the oriented bond sl→target is
 *   - per_bond[sl][n]           if per_bond is non-empty (same shape as
 *                               relative_vectors), used for the anisotropic NN
 *                               bond whose DM matrix varies bond-to-bond, or
 *   - J                         otherwise (a single matrix shared by every bond
 *                               in the shell, e.g. isotropic Heisenberg).
 */
struct CouplingSpec {
    std::string name;
    std::vector<std::vector<ipos_t>> relative_vectors;
    vector3::mat33<double> J;
    std::vector<std::vector<vector3::mat33<double>>> per_bond;
};

struct MC_parameters {
    double T_ref=1.0;
    size_t verbosity = 2;
};

/**
 * @brief Monte Carlo driver for classical spin simulations.
 */
class MC_runner {
    std::vector<CouplingSpec> coupling_specs;
    std::map<std::string, size_t> index;

    vector3::vec3<double> global_field={0,0,0};

    // Single-ion easy-axis anisotropy  H_aniso = -aniso_D * sum_i (S_i . axis)^2.
    // aniso_D > 0 is easy-axis (favours alignment along ±axis); default off.
    double aniso_D = 0.0;
    vector3::vec3d aniso_axis = {0,0,1};

    Lattice* lat;

    std::uniform_int_distribution<size_t> site_dist;
    std::normal_distribution<double> normal_dist;
    std::uniform_real_distribution<double> rand01;

    XoshiroCpp::Xoroshiro128PlusPlus rng;

    vector3::vec3d local_field(const HeisenbergSpin* spin) const;

public:
    MC_parameters settings;

    MC_runner(Lattice &lat_, size_t seed)
        : lat(&lat_),
          site_dist(0, lat_.get_objects<HeisenbergSpin>().size()-1),
          rand01(0,1),
          rng(seed)
    {}

    // Repoints this runner at a different (but congruent) Lattice — O(1), no
    // spin data copied. Used by parallel tempering to swap which replica's
    // configuration a given temperature slot is currently simulating.
    void rebind(Lattice& new_lat);

    double total_energy_per_unit_cell() const;

    // Shell with a single coupling matrix shared by every bond (e.g. Heisenberg).
    void define_coupling(const std::string& name,
            const std::vector<std::vector<ipos_t>>& rel_vecs,
            const vector3::mat33<double>& J);

    // Shell with a per-bond coupling matrix (same shape as rel_vecs). Used for
    // the NN bond whose Dzyaloshinskii–Moriya matrix rotates between bonds.
    void define_coupling(const std::string& name,
            const std::vector<std::vector<ipos_t>>& rel_vecs,
            const std::vector<std::vector<vector3::mat33<double>>>& per_bond);

    void set_global_field(const vector3::vec3<double>& h);
    vector3::vec3d get_global_field() const;

    // Sets the single-ion easy-axis anisotropy -D (S.axis)^2. axis need not be
    // normalised; it is normalised internally.
    void set_anisotropy(double D, const vector3::vec3d& axis);

    void setup_lattice();

    // Single-spin Metropolis update. `sigma` is the standard deviation of the
    // Gaussian kick applied to S before re-normalisation onto the unit sphere.
    size_t local_Metropolis(double T, HeisenbergSpin* spin, double sigma);

    void overrelax_all();
    void overrelax_some(double p);

    // Equilibrium sweep: one Metropolis pass with adaptive proposal width
    // sigma = sqrt(T/T_ref) followed by a microcanonical overrelaxation pass.
    // Decorrelates quickly but is NOT a physical time evolution (overrelaxation
    // instantaneously re-randomises the transverse spin component). Use only for
    // burn-in, annealing, and equilibrium sampling.
    size_t sweep_equilibrate(double T);

    // Pure-Metropolis equilibration sweep: one Metropolis pass over all spins
    // with adaptive proposal width sigma = sqrt(T/T_ref) and NO overrelaxation.
    // Required when single-ion anisotropy is active: overrelaxation reflects
    // spins about the exchange local field and is microcanonical only for the
    // exchange+Zeeman energy, so it would break detailed balance in the
    // presence of the on-site anisotropy term. Pinned spins are skipped.
    size_t sweep_metropolis(double T);

    // Dynamics sweep: N random-sequential single-spin Metropolis attempts
    // (Glauber / Model-A update) with a fixed proposal width `sigma` and no
    // overrelaxation. Random site selection avoids the directional bias of a
    // typewriter sweep. One sweep is then a fixed increment of (overdamped)
    // physical time, suitable for time-dependent measurements such as AC
    // susceptibility. Smaller `sigma` → finer time resolution and a cleaner
    // Langevin correspondence, at the cost of slower evolution per sweep.
    size_t sweep_dynamics(double T, double sigma);
};



// Write a "/geometry" group to an open HDF5 file with lattice statistics
// that are fixed for this disorder realisation:
//
//   n_spins            (scalar)  — number of non-deleted spins
//   n_tetras_by_intact (hsize_t[5]) — n_tetras_by_intact[k] = number of
//                                     tetrahedra with exactly k intact spins
 void write_geometry_group(hid_t file_id, Lattice& sc,
                                 const char* group_name = "geometry");



void save_spin_state(const Lattice& lat, const std::filesystem::path& file_path);


}
