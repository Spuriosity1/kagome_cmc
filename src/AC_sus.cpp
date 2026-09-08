#include <cmath>

#include "MC.hpp"
#include "cli_bits.hpp"
#include "energy_manager.hpp"
#include "format_bits.hpp"

/*
This program performs classical Monte Carlo AC susceptibility sampling for a
spin model on the kagome lattice.

Protocol:
  1. Constructs a supercell and initialises random spins.
  2. Burns in at T_hot.
  3. Anneals from T_hot to T_sample (n_anneal_steps temperature steps).
  4. Applies an oscillating field  B(t) = B_dc + sin(omega * n_sweep * i) * deltaB
     at T_sample, running n_sweep MC sweeps per time step i.
  5. Accumulates M(t) and decomposes into chi'(omega) and chi''(omega).

Output (HDF5):
  /ac_susceptibility/chi_prime[3]   — in-phase susceptibility (per spin)
  /ac_susceptibility/chi_second[3]  — out-of-phase susceptibility (per spin)
  /ac_susceptibility/M_dc[3]        — DC magnetisation (per spin)
  /ac_susceptibility/B[3], deltaB[3], T_sample, omega, n_steps, n_sweep
  /energy                           — annealing energy trace
  /geometry                         — lattice metadata
*/

using namespace std;
using namespace CMC;


int main (int argc, char *argv[]) {
    ///////////////////////////////////////////////////////////////////////////
    /// Setup and CLI
    argparse::ArgumentParser prog("AC_sus");

    /// BOOK-KEEPING
    prog.add_argument("--output_dir", "-o")
        .help("Path to output directory")
        .required();
    prog.add_argument("--seed", "-s")
        .required()
        .help("Seed index into the true-random cache")
        .scan<'i', size_t>();
    prog.add_argument("--save_state")
        .implicit_value(true)
        .default_value(false);
    prog.add_argument("--prefix")
        .default_value("run");

    /// ANNEALING PROTOCOL
    prog.add_argument("--T_hot")
        .help("Temperature to begin annealing from (default: 10)")
        .default_value(10.0)
        .scan<'g', double>();
    prog.add_argument("--T_ref")
        .help("Reference temperature for proposal distribution")
        .default_value(1.0)
        .scan<'g', double>();
    prog.add_argument("--T_sample")
        .help("Sampling temperature (final anneal target)")
        .required()
        .scan<'g', double>();
    prog.add_argument("--n_anneal_steps")
        .help("Number of annealing steps (temperature grid points from T_hot to T_sample)")
        .default_value(static_cast<size_t>(100))
        .scan<'i', size_t>();
    prog.add_argument("--n_burn_in")
        .help("Number of sweeps at T_hot before annealing")
        .default_value(static_cast<size_t>(160))
        .scan<'i', size_t>();

    /// AC SAMPLING
    prog.add_argument("--omega")
        .help("Drive frequency (rad / sweep)")
        .required()
        .scan<'g', double>();
    prog.add_argument("--n_steps")
        .help("Number of AC time steps (each step = n_sweep MC sweeps)")
        .default_value(static_cast<size_t>(1024))
        .scan<'i', size_t>();
    prog.add_argument("--n_sweep", "-w")
        .help("Number of MC sweeps per AC time step")
        .default_value(static_cast<size_t>(16))
        .scan<'i', size_t>();
    prog.add_argument("--deltaB")
        .help("AC field amplitude and direction (3 components)")
        .nargs(3)
        .default_value(std::vector<double>({0., 0., 0.}))
        .scan<'g', double>();

    // PHYSICS
    declare_LJ123(prog);

    try {
        prog.parse_args(argc, argv);
    } catch (const std::exception& err) {
        cerr << err.what() << "\n" << prog;
        std::exit(1);
    }


    ///////////////////////////////////////////////////////////////////////////
    /// Input loading

    filesystem::path outdir(prog.get<std::string>("output_dir"));
    if (!filesystem::exists(outdir))
        throw runtime_error("Cannot open output directory: " + outdir.string());

    const size_t seed         = prog.get<size_t>("--seed");
    const double T_hot        = prog.get<double>("--T_hot");
    const double T_sample     = prog.get<double>("--T_sample");
    const size_t n_anneal     = prog.get<size_t>("--n_anneal_steps");
    const size_t n_burn_in    = prog.get<size_t>("--n_burn_in");
    const double omega        = prog.get<double>("--omega");
    const size_t n_steps      = prog.get<size_t>("--n_steps");
    const size_t n_sweep      = prog.get<size_t>("--n_sweep");

    auto lat = build_lat(prog);
    auto mc  = build_mc_run(prog, lat, seed);

    const auto B_vec     = prog.get<std::vector<double>>("-B");
    const auto dB_vec    = prog.get<std::vector<double>>("--deltaB");
    vector3::vec3d B     = from_std_vector(B_vec);
    vector3::vec3d deltaB = from_std_vector(dB_vec);

    const double deltaB_mag = std::sqrt(deltaB[0]*deltaB[0]
                                       + deltaB[1]*deltaB[1]
                                       + deltaB[2]*deltaB[2]);

    std::stringstream name;
    name << prog.get<std::string>("--prefix") << DELIM
         << name_LJ123(prog)
         << "B="  << B_vec[0] << "," << B_vec[1] << "," << B_vec[2] << DELIM
         << "dB=" << dB_vec[0] << "," << dB_vec[1] << "," << dB_vec[2] << DELIM
         << "om=" << omega << DELIM
         << "seed=" << seed << DELIM
         << "Tc=" << T_sample << DELIM
         << "sw=" << n_sweep << DELIM;

    auto [T_grid, _] = generate_T_profile(T_hot, T_sample, {}, n_anneal);
    energy_manager e_manager;

    auto file_path = outdir / (name.str() + ".out.h5");
    hid_t file_id = H5Fcreate(file_path.string().c_str(), H5F_ACC_TRUNC,
                               H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0)
        throw std::runtime_error("Failed to create HDF5 file: " + file_path.string());

    auto sweep = [&](double T_) -> size_t {
        return mc.sweep_local_Metropolis(T_);
    };

    const size_t n_spins = lat.get_objects<HeisenbergSpin>().size();

    printf("AC_sus: L=%d  n_spins=%zu  omega=%.4e  n_steps=%zu  n_sweep=%zu\n",
           prog.get<int>("L"), n_spins, omega, n_steps, n_sweep);
    printf("        T_hot=%.3g  T_sample=%.3g  n_anneal=%zu  n_burn_in=%zu\n",
           T_hot, T_sample, n_anneal, n_burn_in);
    printf("        deltaB=(%.3g,%.3g,%.3g)  |deltaB|=%.4g\n",
           deltaB[0], deltaB[1], deltaB[2], deltaB_mag);

    // Initialisation
    for (auto& s : lat.get_objects<HeisenbergSpin>()){
        s.S = {1,0,0};
    }

    ///////////////////////////////////////////////////////////////////////////
    /// Stage 0: burn-in at T_hot

    printf("\nBurning in (%zu sweeps at T=%.3g)...\n", n_burn_in, T_hot);
    for (size_t i = 0; i < n_burn_in; i++)
        sweep(T_hot);

    ///////////////////////////////////////////////////////////////////////////
    /// Stage 1: anneal from T_hot to T_sample

    printf("Annealing (%zu steps)...\n", T_grid.size());
    double T = T_hot;
    for (size_t i = 0; i < T_grid.size(); ++i) {
        T = T_grid[i];
        size_t accepted = 0;
        e_manager.new_T(T);
        for (size_t n = 0; n < n_sweep; n++)
            accepted += sweep(T);
        double e = mc.total_energy_per_unit_cell();
        e_manager.sample(e);
        printf("[anneal] %4zu  T=%.3e  E=%.3e  acc=%.1f%%\n",
               i, T, e_manager.curr_E(),
               accepted * 100.0 / n_spins / n_sweep);
    }

    ///////////////////////////////////////////////////////////////////////////
    /// Stage 2: AC susceptibility sampling
    //
    // Field at step i:
    //   B(i) = B_dc + sin(omega * n_sweep * i) * deltaB
    //
    // After N steps the accumulators give:
    //   chi'_alpha  =  (2/N) * Sigma_i  M_alpha(i) * sin(phase_i) / |deltaB|
    //   chi''_alpha = -(2/N) * Sigma_i  M_alpha(i) * cos(phase_i) / |deltaB|
    //   M_dc_alpha  =  (1/N) * Sigma_i  M_alpha(i)
    //
    // The factor 2 arises from the one-sided Fourier decomposition; chi'' has
    // the sign convention that positive chi'' means absorption (energy input).

    printf("\nAC sampling (%zu steps)...\n", n_steps);

    double acc_sin[3] = {0, 0, 0};  // Σ M_α * sin(phase)
    double acc_cos[3] = {0, 0, 0};  // Σ M_α * cos(phase)
    double acc_dc[3]  = {0, 0, 0};  // Σ M_α

    for (size_t i = 0; i < n_steps; i++) {
        const double phase = omega * static_cast<double>(n_sweep) * static_cast<double>(i);
        const double s = std::sin(phase);
        const double c = std::cos(phase);

        mc.set_global_field(B + s * deltaB);

        size_t accepted = 0;
        for (size_t n = 0; n < n_sweep; n++)
            accepted += sweep(T_sample);

        // Measure per-spin magnetisation
        double M[3] = {0, 0, 0};
        for (const auto& spin : lat.get_objects<HeisenbergSpin>()) {
            M[0] += spin.S[0];
            M[1] += spin.S[1];
            M[2] += spin.S[2];
        }
        const double inv_N = 1.0 / static_cast<double>(n_spins);
        M[0] *= inv_N;  M[1] *= inv_N;  M[2] *= inv_N;

        for (int k = 0; k < 3; k++) {
            acc_sin[k] += M[k] * s;
            acc_cos[k] += M[k] * c;
            acc_dc[k]  += M[k];
        }

        if (i % 256 == 0)
            printf("[AC] %5zu/%zu  phase=%.2f  M=(%.3e,%.3e,%.3e)  acc=%.1f%%\n",
                   i, n_steps, phase, M[0], M[1], M[2],
                   accepted * 100.0 / n_spins / n_sweep);
    }

    // Normalise
    const double inv_nsteps = 1.0 / static_cast<double>(n_steps);
    const double ac_norm = (deltaB_mag > 1e-14) ? 2.0 * inv_nsteps / deltaB_mag : 0.0;

    double chi_prime[3], chi_second[3], M_dc[3];
    for (int k = 0; k < 3; k++) {
        chi_prime[k]  =  acc_sin[k] * ac_norm;
        chi_second[k] = -acc_cos[k] * ac_norm;
        M_dc[k]       =  acc_dc[k]  * inv_nsteps;
    }

    printf("\nchi'   = (%.4e, %.4e, %.4e)\n", chi_prime[0],  chi_prime[1],  chi_prime[2]);
    printf("chi''  = (%.4e, %.4e, %.4e)\n",   chi_second[0], chi_second[1], chi_second[2]);
    printf("M_dc   = (%.4e, %.4e, %.4e)\n",   M_dc[0],       M_dc[1],       M_dc[2]);

    ///////////////////////////////////////////////////////////////////////////
    /// Write results

    {
        hid_t grp = H5Gcreate2(file_id, "/ac_susceptibility",
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (grp < 0)
            throw std::runtime_error("Failed to create /ac_susceptibility group");

        const hsize_t dim3 = 3;
        auto write_vec3 = [&](const char* dname, const double v[3]) {
            hid_t sp = H5Screate_simple(1, &dim3, nullptr);
            hid_t ds = H5Dcreate2(grp, dname, H5T_NATIVE_DOUBLE, sp,
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            if (ds < 0) { H5Sclose(sp); throw std::runtime_error(
                std::string("Failed to create dataset ") + dname); }
            H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v);
            H5Dclose(ds);
            H5Sclose(sp);
        };

        auto write_scalar_d = [&](const char* dname, double v) {
            hid_t sp = H5Screate(H5S_SCALAR);
            hid_t ds = H5Dcreate2(grp, dname, H5T_NATIVE_DOUBLE, sp,
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            if (ds < 0) { H5Sclose(sp); throw std::runtime_error(
                std::string("Failed to create dataset ") + dname); }
            H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v);
            H5Dclose(ds);
            H5Sclose(sp);
        };

        auto write_scalar_u = [&](const char* dname, size_t v) {
            hid_t sp = H5Screate(H5S_SCALAR);
            hid_t ds = H5Dcreate2(grp, dname, H5T_NATIVE_ULLONG, sp,
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            if (ds < 0) { H5Sclose(sp); throw std::runtime_error(
                std::string("Failed to create dataset ") + dname); }
            H5Dwrite(ds, H5T_NATIVE_ULLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v);
            H5Dclose(ds);
            H5Sclose(sp);
        };

        const double B_arr[3]  = {B[0],     B[1],     B[2]};
        const double dB_arr[3] = {deltaB[0], deltaB[1], deltaB[2]};

        write_vec3("chi_prime",  chi_prime);
        write_vec3("chi_second", chi_second);
        write_vec3("M_dc",       M_dc);
        write_vec3("B",          B_arr);
        write_vec3("deltaB",     dB_arr);
        write_scalar_d("T_sample",    T_sample);
        write_scalar_d("omega",       omega);
        write_scalar_d("deltaB_mag",  deltaB_mag);
        write_scalar_u("n_steps",     n_steps);
        write_scalar_u("n_sweep",     n_sweep);

        H5Gclose(grp);
    }

    e_manager.write_group(file_id, "/energy");
    write_geometry_group(file_id, lat);
    H5Fclose(file_id);
    std::cout << "Saved to " << file_path << "\n";

    if (prog.get<bool>("--save_state")) {
        auto f = outdir / (name.str() + ".spins.h5");
        printf("Saving spin state to %s\n", f.string().c_str());
        save_spin_state(lat, f);
    }

    return 0;
}
