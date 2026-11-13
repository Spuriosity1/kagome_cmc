#pragma once
#include "format_bits.hpp"
#include "vec3.hpp"
#include <algorithm>
#include <array>
#include <argparse/argparse.hpp>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <ostream>

#include "vec3.hpp"
#include "random_seeds.hpp"
#include "geometry/kagome.hpp"



using mat33d=vector3::mat33<double>;

namespace coupling {
static const mat33d Heis = mat33d::from_cols(
        {1,0,0},{0,1,0},{0,0,1});

static const mat33d Isin = mat33d::from_cols(
        {0,0,0},{0,0,0},{0,0,1});
};


inline auto declare_LJ123(argparse::ArgumentParser& prog){

    /// PHYSICAL
    prog.add_argument("L")
        .help("Linear dimension of cubic supercell (e.g. L=2 has 2^3 *4 = 32 primitive cells)")
        .scan<'i', int>(); // cubic unit cell hardcoded
    prog.add_argument("--J1")
        .help("Nearest-neighbour Heisenberg coupling strength")
        .default_value(0.)
        .scan<'g', double>();
    prog.add_argument("--J2")
        .help("Second-nearest-neighbour Heisenberg coupling strength")
        .default_value(0.)
        .scan<'g', double>();
    prog.add_argument("--J3")
        .help("Third-nearest-neighbour Heisenberg coupling strength across bonds")
        .default_value(0.)
        .scan<'g', double>();
    prog.add_argument("--Jd")
        .help("Third-nearest-neighbour Heisenberg coupling strength across hexagons")
        .default_value(0.)
        .scan<'g', double>();
    prog.add_argument("--Jc1")
        .help("Out-of-plane Heisenberg coupling strength")
        .default_value(0.)
        .scan<'g', double>();
    prog.add_argument("--external_field", "-B")
        .help("Global magnetic field")
        .nargs(3)
        .default_value(std::vector<double>({0.,0.,0.}))
        .scan<'g', double>();
}

// Declares the nearest-neighbour Dzyaloshinskii-Moriya options. With these set
// the NN (J1) bonds carry the symmetry-adapted DM of the kagome layer: a
// Cartesian DM vector that rotates by the C3 symmetry between bonds (see
// kagome::dm_J1). It is parametrised by an in-plane strength D1p and an
// out-of-plane strength D1z, which the CLI flags map onto as
//     D1p = -D,   D1z = d*D,
// i.e. --D sets the overall DM scale and --d the out-of-plane / in-plane ratio.
//
// NOTE: the domain_wall protocol's single-ion anisotropy uses a separate flag,
// --D_sia (see declare_cuboid_and_aniso), so there is no collision with --D;
// declare_DM is used by anneal / AC_sus / ltgs and by domain_wall, which shares
// anneal's J1+DM model (only the supercell shape differs).
inline void declare_DM(argparse::ArgumentParser& prog){
    prog.add_argument("--D")
        .help("Nearest-neighbour Dzyaloshinskii-Moriya strength (in-plane "
              "component D1p = -D). Default 0.")
        .default_value(0.)
        .scan<'g', double>();
    prog.add_argument("--d")
        .help("Dimensionless DM anisotropy ratio; out-of-plane component "
              "D1z = d*D. Default 0.")
        .default_value(0.)
        .scan<'g', double>();
}

// Cartesian DM matrix: S_i . dm_matrix(D) . S_j = D . (S_i × S_j).
inline mat33d dm_matrix(const vector3::vec3d& D){
    return mat33d::from_cols({0., -D[2], D[1]},
                             {D[2], 0., -D[0]},
                             {-D[1], D[0], 0.});
}

// Per-bond NN coupling matrices, parallel to kagome::vecs["J1"]: isotropic J1
// plus the symmetry-adapted DM of kagome::dm_J1 (CLI mapping D1p=-D, D1z=d*D).
// The oriented bond sl→target carries H_bond = J1 S_sl·S_t + D·(S_sl×S_t).
inline std::vector<std::vector<mat33d>>
nn_coupling_per_bond(double J1, double D, double d){
    const double D1p = -D, D1z = d * D;
    const auto dm = kagome::dm_J1(D1p, D1z);
    std::vector<std::vector<mat33d>> out(dm.size());
    for (size_t sl = 0; sl < dm.size(); sl++)
        for (const auto& Dv : dm[sl])
            out[sl].push_back(J1 * coupling::Heis + dm_matrix(Dv));
    return out;
}

// Filename fragment for the DM parameters (only when explicitly set).
inline std::string name_DM(const argparse::ArgumentParser& prog){
    std::stringstream s;
    if (prog.is_used("--D")) s << "D=" << prog.get<double>("--D") << DELIM;
    if (prog.is_used("--d")) s << "d=" << prog.get<double>("--d") << DELIM;
    return s.str();
}

// Declares the single-ion easy-axis anisotropy options for the anneal / AC
// protocols:  H += -D_sia * sum_i (S_i . axis)^2.  Default off (D_sia=0) so
// existing runs are unchanged. These share the --D_sia / --easy_axis names used
// by the domain_wall protocol (declare_cuboid_and_aniso) but, since each binary
// builds its own parser, there is no collision; they are kept separate because
// anneal does not want domain_wall's cuboid / pinning flags. NOTE: the on-site
// term breaks the
// microcanonical overrelaxation step's detailed balance, so anneal must use the
// pure-Metropolis sweep (sweep_metropolis), which it already does.
inline void declare_single_ion_aniso(argparse::ArgumentParser& prog){
    prog.add_argument("--D_sia")
        .help("Single-ion easy-axis anisotropy strength (H += -D_sia (S.axis)^2). "
              "D_sia>0 is easy-axis. Default 0 (off).")
        .default_value(0.).scan<'g', double>();
    prog.add_argument("--easy_axis")
        .help("Easy axis direction (3 components, need not be normalised). "
              "Default (0,0,1).")
        .nargs(3)
        .default_value(std::vector<double>({0.,0.,1.}))
        .scan<'g', double>();
}

// Filename fragment for the single-ion anisotropy (only when explicitly set).
inline std::string name_aniso(const argparse::ArgumentParser& prog){
    std::stringstream s;
    if (prog.is_used("--D_sia")) s << "Dsia=" << prog.get<double>("--D_sia") << DELIM;
    return s.str();
}

inline uint64_t int_from_hex_str(const std::string& seed_s){
    uint64_t seed; // ugly hack for loading seed as hex
    std::stringstream ss;
    ss << std::hex << seed_s;
    ss >> seed; 
    return seed;
}

inline auto build_lat(const argparse::ArgumentParser& prog){
    int L = prog.get<int>("L");
    // Orthogonal L×L×L supercell (in primitive-cell units). The interlayer ABC
    // shift is carried by the primitive cell's a3, so no supercell shear is
    // needed — this matches the notebook's simple N1×N2×N3 fractional supercell.
    auto supercell_spec = imat33_t::from_cols({2*L,0,0},{0,2*L,0},{0,0,L});
    auto cell_spec = KagomeLayer();
    return build_supercell(cell_spec, supercell_spec);
}

// Declares the extra options used by the domain_wall protocol: independent
// cuboid dimensions (each defaulting to the cubic L when unset), the by-hand
// single-ion easy-axis anisotropy strength/axis, and the pinned-slab thickness.
inline void declare_cuboid_and_aniso(argparse::ArgumentParser& prog){
    prog.add_argument("--Lx")
        .help("Linear dimension along x (wall normal). Defaults to L.")
        .default_value(-1).scan<'i', int>();
    prog.add_argument("--Ly")
        .help("Linear dimension along y. Defaults to L.")
        .default_value(-1).scan<'i', int>();
    prog.add_argument("--Lz")
        .help("Linear dimension along z (stacking axis). Defaults to L.")
        .default_value(-1).scan<'i', int>();
    prog.add_argument("--D_sia")
        .help("Single-ion easy-axis anisotropy strength (H += -D_sia (S.axis)^2). "
              "D_sia>0 is easy-axis. Default 1.36 (physical Co3Sn2S2 value).")
        .default_value(1.36).scan<'g', double>();
    prog.add_argument("--easy_axis")
        .help("Easy axis direction (3 components, need not be normalised)")
        .nargs(3)
        .default_value(std::vector<double>({0.,0.,1.}))
        .scan<'g', double>();
    prog.add_argument("--pin_thickness")
        .help("Thickness (in integer x units) of each pinned end slab. Must be "
              ">= the max bond range in x (=2) for clean open-BC behaviour.")
        .default_value(2).scan<'i', int>();
}

// Builds an ORTHOGONAL cuboid supercell {2Lx,0,0},{0,2Ly,0},{0,0,Lz}. Unlike
// build_lat (which shears the second vector), this keeps transverse wrapping
// from shifting x, so an x-resolved domain-wall profile is well defined.
inline auto build_cuboid_lat(const argparse::ArgumentParser& prog){
    int L  = prog.get<int>("L");
    int Lx = prog.get<int>("--Lx"); if (Lx < 0) Lx = L;
    int Ly = prog.get<int>("--Ly"); if (Ly < 0) Ly = L;
    int Lz = prog.get<int>("--Lz"); if (Lz < 0) Lz = L;
    auto supercell_spec = imat33_t::from_cols({2*Lx,0,0},{0,2*Ly,0},{0,0,Lz});
    auto cell_spec = KagomeLayer();
    return build_supercell(cell_spec, supercell_spec);
}

// Convenience: resolve the (possibly-defaulted) cuboid dimensions.
inline std::array<int,3> get_cuboid_dims(const argparse::ArgumentParser& prog){
    int L  = prog.get<int>("L");
    int Lx = prog.get<int>("--Lx"); if (Lx < 0) Lx = L;
    int Ly = prog.get<int>("--Ly"); if (Ly < 0) Ly = L;
    int Lz = prog.get<int>("--Lz"); if (Lz < 0) Lz = L;
    return {Lx, Ly, Lz};
}

inline uint64_t hash_true_random(uint64_t seed_id){
    if (seed_id > true_random_cache.size()){
        throw std::runtime_error("Cannot read: there are only "+std::to_string(true_random_cache.size())+"seeds");
    }
    return true_random_cache[seed_id];
}

inline std::map<std::string, double> get_J(const argparse::ArgumentParser& prog){
    std::map<std::string, double> J;

    for (const auto& [cn,_] : kagome::vecs){
        if (prog.is_used(cn)){
            J[cn] = prog.get<double>(cn);
        }
    }

    return J;
}


// D_dm / d_dm add the symmetry-adapted Dzyaloshinskii-Moriya term to the NN
// (J1) bonds (default off). Each oriented NN bond carries its own DM matrix
// (nn_coupling_per_bond); setup_lattice gives the source spin that matrix and
// the target spin its transpose, so the (non-symmetric) bond is felt correctly
// from both ends.
inline auto build_mc_run(const argparse::ArgumentParser& prog, CMC::Lattice& lat,
        uint64_t seed, double D_dm = 0.0, double d_dm = 0.0){

    vector3::vec3d global_field;
    {
        auto B_tmp = prog.get<std::vector<double>>("-B");
        for (int i=0; i<3; i++){ global_field[i]= B_tmp[i]; }
    }


    CMC::MC_runner mc(lat, hash_true_random(seed));

    mc.set_global_field(global_field);

    auto J = get_J(prog);
    const double J1val = J.count("J1") ? J.at("J1") : 0.0;
    const bool dm_active = (std::abs(D_dm) > 1e-12 || std::abs(d_dm * D_dm) > 1e-12);

    for (const auto& [cn, v] : J ){
        if (cn == "J1") continue; // NN coupling handled below (may carry DM)
        if (abs(v) > 1e-10) mc.define_coupling(cn, kagome::vecs.at(cn), v*coupling::Heis);
    }
    // Nearest-neighbour coupling: isotropic J1 plus the (optional) symmetry DM.
    if (dm_active || abs(J1val) > 1e-10){
        mc.define_coupling("J1", kagome::vecs.at("J1"),
                           nn_coupling_per_bond(J1val, D_dm, d_dm));
    }

    mc.settings.T_ref = prog.get<double>("--T_ref");
    mc.setup_lattice();

    return mc;
}

inline auto name_LJ123(const argparse::ArgumentParser& prog){
    int L = prog.get<int>("L");
    auto J = get_J(prog);

    std::stringstream name;
    name << "L="<<L<<DELIM;
    for (auto& [k,v] : J){
        name<<k<<"="<<v<<DELIM;
    }
    return name.str();
}

// Generates a log-spaced array from T_hot to min(T_sample), making sure to
// include all of the T_sample steps exactly in sorted order.
// Returns also the 
inline std::pair<std::vector<double>, std::set<size_t>> generate_T_profile(
        double T_hot, double T_cold,
        std::vector<double> T_sample, size_t total_steps){
    std::vector<double> T_sweep, T;
    const double factor = pow(T_cold / T_hot, 1./total_steps);
    double t = T_hot/factor;
    for (size_t i=0; i<total_steps; i++){
        t *= factor;
        T_sweep.push_back(t);
    }

    // sort both descending
    std::sort(T_sweep.begin(),T_sweep.end(), [](double a, double b){return a>b;}); 
    std::sort(T_sample.begin(),T_sample.end(), [](double a, double b){return a>b;}); 
    
    size_t i_samp=0;
    size_t i_sw=0;

    std::set<size_t> idx;
    while(i_sw < T_sweep.size() && i_samp < T_sample.size()){
        if (T_sweep[i_sw] > T_sample[i_samp]){
            T.push_back(T_sweep[i_sw]); i_sw++;
        } else {
            idx.insert(T.size());
            T.push_back(T_sample[i_samp]); i_samp++;
        }
    }
    
    // exhaust any tails
    for (; i_samp < T_sample.size(); i_samp++){
        idx.insert(T.size());
        T.push_back(T_sample[i_samp]);
    }
    for (; i_sw < T_sweep.size(); i_sw++)
        T.push_back(T_sweep[i_sw]);


    return {T, idx};

}


