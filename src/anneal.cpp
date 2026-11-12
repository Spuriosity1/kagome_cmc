#include <algorithm>

#include "MC.hpp"
#include "cli_bits.hpp"
#include "ssf_manager.hpp"
#include "energy_manager.hpp"
#include "order_param_manager.hpp"
#include "geometry/kagome.hpp"
#include "format_bits.hpp"

/*
This program performs classical Monte Carlo simulated annealing for a spin
model on the kagome lattice, using local Metropolis updates. The simulation:

1. Constructs a supercell of the kagome lattice.
2. Defines exchange couplings up to third neighbors (plus interlayer Jc1).
3. Performs a burn-in at high temperature.
4. Anneals the system from a hot temperature to a target cold temperature.
5. Samples static spin–spin correlation functions at the sample temperatures.
6. Optionally saves the final spin configuration.

The output consists of:

a. Static structure factor data, written to the /ssf group.
b. Optionally, the final spin state (.spins.h5)

 */

using namespace std;
using namespace CMC;


int main (int argc, char *argv[]) {
    ///////////////////////////////////////////////////////////////////////////
    /// Setup and CLI
    argparse::ArgumentParser prog("anneal");

    /// BOOK-KEEPING
    prog.add_argument("--output_dir", "-o")
        .help("Path to output")
        .required();

    prog.add_argument("--seed", "-s")
        .required()
        .help("Seed index to seed the RNG")
        .scan<'i', size_t>();

    prog.add_argument("--save_state")
        .implicit_value(true)
        .default_value(false);

    prog.add_argument("--init_spiral")
        .help("Pre-initialise spins to a spiral with the wavevector given by --Q (requires --Q)")
        .implicit_value(true)
        .default_value(false);


    /// ANNEALING PROTOCOL
    prog.add_argument("--T_hot")
        .help("Temperature  to begin annealing from")
        .scan<'g',double>();
    prog.add_argument("--T_cold")
        .help("Minimum temperature (defaults ti smallest T_sample)")
        .scan<'g',double>();
    prog.add_argument("--T_ref")
        .help("Reference temperature for proposal distribution")
        .default_value(1.0)
        .scan<'g',double>();
    prog.add_argument("--T_sample")
        .help("Final temperature")
        .nargs(argparse::nargs_pattern::at_least_one)
        .required()
        .scan<'g',double>();
    prog.add_argument("--n_steps")
        .help("Number of annealing steps")
        .default_value(static_cast<size_t>(100))
        .scan<'i', size_t>();
    prog.add_argument("--n_sweep", "-w")
        .help("Number of sweeps to run per temperatur step")
        .default_value(static_cast<size_t>(16))
        .scan<'i', size_t>();
    prog.add_argument("--n_burn_in")
        .help("Number of sweeps to run at T_hot before we start annealing")
        .default_value(static_cast<size_t>(16))
        .scan<'i', size_t>();
    prog.add_argument("--n_sample")
        .default_value(static_cast<size_t>(64))
        .help("Number of sweeps to run at T_cold while collecting statistics")
        .scan<'i', size_t>();
    prog.add_argument("--n_equil")
        .help("Number of sweeps to equilibrate at each new temperature before "
              "collecting statistics (per-T burn-in). Too small and the "
              "specific heat is inflated by the cooling transient.")
        .default_value(static_cast<size_t>(256))
        .scan<'i', size_t>();
    prog.add_argument("--prefix")
        .default_value("run");

    declare_LJ123(prog);
    declare_DM(prog);
    declare_single_ion_aniso(prog);

    try {
        prog.parse_args(argc, argv);
    } catch (const std::exception& err){
        cerr << err.what() << endl;
        cerr << prog;
        std::exit(1);
    }


    ///////////////////////////////////////////////////////////////////////////
    /// Input loading and validation
    

    /// Ensuring directories exist AHEAD of time (avoids heartbreak)
    std::string outdir_s = prog.get<std::string>("output_dir");
    filesystem::path outdir(outdir_s);
    if (! filesystem::exists(outdir) ){
        throw runtime_error("Cannot open outdir");
    }
    
    size_t seed = prog.get<size_t>("--seed");

    const double T_hot = prog.is_used("--T_hot") ? 
        prog.get<double>("--T_hot") : 10;
    const std::vector<double> T_sample = prog.get<std::vector<double>>("--T_sample");
    double T_cold =prog.is_used("--T_cold") ? 
        prog.get<double>("--T_cold") : *std::min_element(T_sample.begin(), T_sample.end());

    const size_t n_steps = prog.get<size_t>("--n_steps");
    const size_t n_sweep = prog.get<size_t>("--n_sweep");
    const size_t n_burn_in = prog.get<size_t>("--n_burn_in");
    const size_t n_sample = prog.get<size_t>("--n_sample");
    const size_t n_equil = prog.get<size_t>("--n_equil");

    double T = T_hot;

    auto lat = build_lat(prog);
    auto mc = build_mc_run(prog, lat, seed,
                           prog.get<double>("--D"), prog.get<double>("--d"));

    // Single-ion easy-axis anisotropy: H += -D_sia (S.axis)^2. Off (D_sia=0) by
    // default. anneal sweeps with sweep_metropolis, which keeps detailed balance
    // under this on-site term.
    const double D_sia = prog.get<double>("--D_sia");
    const vector3::vec3d easy_axis =
        from_std_vector(prog.get<std::vector<double>>("--easy_axis"));
    mc.set_anisotropy(D_sia, easy_axis);

    auto B = prog.get<std::vector<double>>("-B");
    if (B.size() < 3){
        throw std::runtime_error("--B requires 3 arguments");
    }
    // Parameter specification complete. Set the name...
    std::stringstream name; // accumulates hashed options
    name << prog.get<std::string>("--prefix")<<DELIM<<name_LJ123(prog)<<name_DM(prog)<<name_aniso(prog)<<
        "B="<<B[0]<<","<<B[1]<<","<<B[2]<<DELIM<<
        "seed="<<seed<<DELIM<<
        "Tc="<<T_cold<<DELIM<<
        "sw="<<n_sweep<<DELIM;

    auto [T_grid, S_idx] = generate_T_profile(T_hot, T_cold, T_sample, n_steps);
    energy_manager e_manager;
    order_manager  o_manager(lat);

    // Create the output file before the sampling loop so ssf_manager can
    // pre-allocate static_corr on disk and write each temperature's data
    // immediately after sampling (streaming mode).  File creation here also
    // means a partial result is recoverable if the job is killed mid-run.
    auto file_path = outdir/( name.str() + ".out.h5");
    hid_t file_id = H5Fcreate(file_path.string().c_str(), H5F_ACC_TRUNC,
                               H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0)
        throw std::runtime_error("Failed to create HDF5 file: " + file_path.string());

    ssf_manager ssfm(lat, {"xx", "yy", "zz"}, file_id, "/ssf", T_sample, true);

    auto sweep = [&](double T_) -> size_t {
        return mc.sweep_metropolis(T_);
    };

    printf("Done. Begin anneal...\n");

    printf("Burning in (%zu sweeps)...\n", n_burn_in);
    for (size_t i=0; i<n_burn_in; i++){
        sweep(T_hot);
    }

    for (size_t i=0; i<T_grid.size(); ++i){
        T = T_grid[i];
        e_manager.new_T(T);
        o_manager.new_T(T);

        // Equilibrate at this temperature BEFORE measuring. The configuration
        // arrives still relaxed to the previous (hotter) step; sampling while
        // it relaxes injects the cooling transient into Var(E) and inflates the
        // specific heat. This per-T burn-in is what lets each fixed-T chain
        // converge before we read off statistics.
        for (size_t n=0; n<n_equil; n++){
            sweep(T);
        }

        size_t accepted_Cv=0;
        for (size_t s=0; s<n_sample; s++){
            // Decorrelate first, then measure — never sample the stale
            // inherited configuration.
            for (size_t n=0; n<n_sweep; n++){
                accepted_Cv += sweep(T);
            }
            e_manager.sample(mc.total_energy_per_unit_cell());
            o_manager.sample();
        }

        if (S_idx.contains(i)){
            printf("Sampling at T=%lf (%zu sweeps)...\n", T, n_sample);
            ssfm.new_T(T);
            for (size_t i=0; i<n_sample; i++){
                size_t accepted_ssf=0;
                for (size_t n=0; n<n_sweep; n++){
                    accepted_ssf += sweep(T);
                }
                printf("[ssf] Iter %4zu T=%.3e Acceptance rate: %.2f%%\n",
                        i, T, accepted_ssf*100.0/lat.get_objects<HeisenbergSpin>().size()/n_sweep);
                ssfm.sample();
            }
            ssfm.flush();
        }

        double E = e_manager.curr_E();
        printf("[anneal] Iter %4zu T=%.3e E=%3e Acceptance rate: %.2f%%\n", 
                i, T, E, accepted_Cv*100.0/lat.get_objects<HeisenbergSpin>().size()/(
                    n_sweep * n_sample) );

    }

    ssfm.write_group(file_id, "/ssf");
    e_manager.write_group(file_id, "/energy");
    o_manager.write_group(file_id, "/order");
    write_geometry_group(file_id, lat);
    H5Fclose(file_id);
    std::cout<<"Saved to \n"<< file_path<<std::endl;
    
    if (prog.get<bool>("--save_state")){
        auto f = outdir /( name.str() + ".spins.h5" );
        printf("Saving spin state to %s\n", f.string().c_str());
        save_spin_state(lat, f);
    }

    return 0;
}
