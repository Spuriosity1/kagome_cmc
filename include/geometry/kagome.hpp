#pragma once

#include "MC.hpp"


/**
 * Factory function: builds a UnitCellSpecifier for the kagome spin lattice.
 *
 */
inline UnitCellSpecifier<CMC::HeisenbergSpin> KagomeLayer() {
    UnitCellSpecifier<CMC::HeisenbergSpin> spec(
        imat33_t::from_cols({2,0,0}, {0,2,0}, {0,0,1}));

	static const std::vector<ipos_t> sl_pos = {
		{0,0,0}, {1,0,0}, {1,1,0}
	};

	for (int p = 0; p < 3; p++) {
		CMC::HeisenbergSpin spin;
		spin.ipos = sl_pos[p];
		spin.pyro_sl = p;
		spec.add(std::move(spin));
	}
    
    return spec;
}

namespace kagome {

    using idx3_t = ipos_t;
    using v_sl = std::vector<std::vector<idx3_t>>;

    const std::map<std::string, v_sl> vecs = {
    {"J1", {
        {{1,0,0},{-1,0,0}}, // SL 0 @ (0,0,0)
        {{0,1,0},{0,-1,0}}, // SL 1 @ (1,0,0)
        {{1,1,0},{-1,-1,0}}  // SL 2 @ (1,1,0)
    }},

    {"J2", {
        {{1,2,0},{-1,1,0}},
        {{1,2,0},{2,1,0} },
        {{-1,1,0},{2,1,0}}
    }},

    {"J3", {
        {{2,0,0},{2,2,0}},
        {{2,0,0},{0,2,0}},
        {{0,2,0},{2,2,0}}
    }},

    {"Jd", {
        {{0,2,0}},
        {{2,2,0}},
        {{2,0,0}}
    }},

    // Interlayer coupling, first out-of-plane shell. Co3Sn2S2 is rhombohedral
    // (R-3m): successive kagome layers are ABC-stacked, each shifted in-plane by
    // the rhombohedral vector +-(2/3 a1 + 1/3 a2). An interlayer bond therefore
    // connects DIFFERENT sublattices (a naive AA c-axis bond {0,0,1} would wrongly
    // link a site to its own sublattice). Each site has two nearest neighbours in
    // the layer above (listed here) and two below (the reverse bonds, supplied
    // automatically). Only the bond TOPOLOGY is physical: the integer embedding
    // keeps every layer at the same in-plane coordinates (a shear of the true ABC
    // structure), preserving the connectivity graph but not the metric.
    {"Jc1", {
        {{1,0,1},{1,1,1}}, // SL 0 @ (0,0,0) -> SL1, SL2 above
        {{1,0,1},{2,1,1}}, // SL 1 @ (1,0,0) -> SL0, SL2 above
        {{1,1,1},{2,1,1}}  // SL 2 @ (1,1,0) -> SL0, SL1 above
    }}
    };

} // end namespace kagome
