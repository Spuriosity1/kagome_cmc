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

    {"Jc1", {
        {{0,0,1}},
        {{0,0,1}},
        {{0,0,1}}
    }}
    };

} // end namespace kagome
