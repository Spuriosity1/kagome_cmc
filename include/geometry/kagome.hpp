#pragma once

#include "MC.hpp"

#include <cmath>


// ---------------------------------------------------------------------------
// Kagome spin lattice.
//
// This geometry is a direct port of the hand-checked Luttinger–Tisza notebook
// (nb/LuttingerTisza.ipynb). The notebook works in fractional coordinates with
// lattice vectors a1=(2,0,0), a2=(1,√3,0), a3=(0,2/√3,2) and sublattices at
// fractional (0,0,0), (½,0,0), (0,½,0). We embed it on the integer lattice by
//     integer = (2·f_x, 2·f_y, f_z),
// so the primitive cell is diag(2,2,1) and the three sublattices sit at integer
// (0,0,0), (1,0,0), (0,1,0). The embedding is an exact relabelling of the
// notebook lattice: the bond topology, the Fourier spectrum and the (Cartesian)
// DM vectors are identical to the notebook by construction.
//
// Bond convention: every undirected bond is listed exactly ONCE (from one
// endpoint). The reverse bond is supplied automatically — by the MC lattice
// builder (setup_lattice adds the Hermitian partner to the other endpoint) and
// by the Luttinger–Tisza code (build_bonds emits β→α alongside α→β). Do NOT
// list a bond from both ends or it will be double counted.
// ---------------------------------------------------------------------------

namespace kagome {

    using idx3_t = ipos_t;
    using v_sl   = std::vector<std::vector<idx3_t>>;

    // Integer sublattice positions in the diag(2,2,1) primitive cell.
    inline const std::vector<ipos_t> sl_pos = {
        {0,0,0}, {1,0,0}, {0,1,0}
    };

    // Per-sublattice bond displacement tables (integer coordinates).
    //   J1  : nearest neighbour (triangle edges), distance a.
    //   J2  : second neighbour, distance √3·a.
    //   J3  : third neighbour "across a bond" (collinear through a NN site),
    //         distance 2a; connects a sublattice to itself.
    //   Jd  : third neighbour "across a hexagon" (through the hexagon void),
    //         distance 2a; connects a sublattice to itself.
    //   Jc1 : first interlayer shell (one layer up). Co3Sn2S2 is rhombohedral
    //         (R-3m) ABC kagome stacking; the in-plane offset between layers is
    //         carried by a3's in-plane component, so an interlayer bond connects
    //         DIFFERENT sublattices. Each site has two neighbours in the layer
    //         above (listed) and two below (the automatic reverse bonds).
    //
    // J1 and Jc1 reproduce the notebook's bond lists exactly (the J1 ordering is
    // also the ordering of the DM table dm_J1 below). J2/J3/Jd are the further
    // Heisenberg shells, derived geometrically in this decoration.
    const std::map<std::string, v_sl> vecs = {
    {"J1", {
        {{1,0,0},{0,1,0}},    // SL0
        {{1,0,0},{1,-1,0}},   // SL1
        {{0,1,0},{1,-1,0}}    // SL2
    }},

    {"J2", {
        {{1,-2,0},{2,-1,0}},  // SL0
        {{1,-2,0},{1,1,0}},   // SL1
        {{1,1,0},{2,-1,0}}    // SL2
    }},

    {"J3", {
        {{0,2,0},{2,0,0}},    // SL0
        {{2,-2,0},{2,0,0}},   // SL1
        {{0,2,0},{2,-2,0}}    // SL2
    }},

    {"Jd", {
        {{2,-2,0}},           // SL0
        {{0,2,0}},            // SL1
        {{2,0,0}}             // SL2
    }},

    {"Jc1", {
        {{0,-1,1},{1,-2,1}},  // SL0 -> SL2, SL1 above
        {{1,-1,1},{1,-2,1}},  // SL1 -> SL2, SL0 above
        {{1,-1,1},{0,-1,1}}   // SL2 -> SL1, SL0 above
    }}
    };

    // Nearest-neighbour Dzyaloshinskii–Moriya vectors, in Cartesian spin
    // coordinates, one per oriented J1 bond (parallel to vecs.at("J1")). This is
    // the notebook's symmetry_dm("J1") with in-plane strength D1p and
    // out-of-plane strength D1z; the C3 symmetry rotates the in-plane part by
    // 120° between bonds (s = √3/2). The oriented bond sl -> target carries
    //     H_bond = J1 (S_sl · S_t) + D · (S_sl × S_t).
    inline std::vector<std::vector<vector3::vec3d>> dm_J1(double D1p, double D1z) {
        const double s = std::sqrt(3.0) / 2.0;
        return {
            {{0.0,    D1p,       D1z}, {-s*D1p,  0.5*D1p,  -D1z}},  // SL0
            {{0.0,   -D1p,      -D1z}, {-s*D1p, -0.5*D1p,   D1z}},  // SL1
            {{s*D1p, -0.5*D1p,   D1z}, { s*D1p,  0.5*D1p,  -D1z}}   // SL2
        };
    }

} // end namespace kagome


/**
 * Factory function: builds a UnitCellSpecifier for the kagome spin lattice.
 */
inline UnitCellSpecifier<CMC::HeisenbergSpin> KagomeLayer() {
    UnitCellSpecifier<CMC::HeisenbergSpin> spec(
        imat33_t::from_cols({2,0,0}, {0,2,0}, {0,0,1}));

    for (int p = 0; p < 3; p++) {
        CMC::HeisenbergSpin spin;
        spin.ipos    = kagome::sl_pos[p];
        spin.pyro_sl = p;
        spec.add(std::move(spin));
    }

    return spec;
}
