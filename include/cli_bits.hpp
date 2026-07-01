#pragma once
#include "format_bits.hpp"
#include "vec3.hpp"
#include <algorithm>
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

inline uint64_t int_from_hex_str(const std::string& seed_s){
    uint64_t seed; // ugly hack for loading seed as hex
    std::stringstream ss;
    ss << std::hex << seed_s;
    ss >> seed; 
    return seed;
}

inline auto build_lat(const argparse::ArgumentParser& prog){
    int L = prog.get<int>("L");
    auto supercell_spec = imat33_t::from_cols({2*L,0,0},{L,2*L,0},{0,0,L});
    auto cell_spec = KagomeLayer();
    return build_supercell(cell_spec, supercell_spec);
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


inline auto build_mc_run(const argparse::ArgumentParser& prog, CMC::Lattice& lat, uint64_t seed){

    vector3::vec3d global_field;
    {
        auto B_tmp = prog.get<std::vector<double>>("-B");
        for (int i=0; i<3; i++){ global_field[i]= B_tmp[i]; }
    }


    CMC::MC_runner mc(lat, hash_true_random(seed));

    mc.set_global_field(global_field);

    auto J = get_J(prog);

    for (const auto& [cn, v] : J ){
        if (abs(v) > 1e-10) mc.define_coupling(cn, kagome::vecs.at(cn), v*coupling::Heis);
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


