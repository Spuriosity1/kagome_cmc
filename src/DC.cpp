#include <array>
#include <cmath>
#include <vector>

#include "MC.hpp"
#include "cli_bits.hpp"
#include "energy_manager.hpp"
#include "format_bits.hpp"

/*
This program performs classical Monte Carlo in a STATIC ("DC") field to study a
magnetic domain wall on the stacked-kagome lattice (parameterised for the kagome
ferromagnet Co3Sn2S2).

Geometry:
  An ORTHOGONAL cuboid supercell elongated along x (the wall normal). The easy
  axis is c (z). The two end slabs in x (thickness --pin_thickness integer
  units) are pinned to +z and -z respectively, imposing Dirichlet boundary
  conditions that stabilise a single domain wall; y and z remain periodic.
  Because the maximum bond range in x is 2 integer units, a pinned slab of
  thickness >= 2 makes the bulk behave as open-BC: wrap-around bonds across the
  x-seam only ever connect frozen spins.

Hamiltonian:
  H = sum_{<ij>_n} J_n S_i . S_j  -  B . sum_i S_i  -  D_sia sum_i (S_i . axis)^2
  with the single-ion easy-axis anisotropy D_sia added by hand (--D_sia) so the discrete
  wall can be pushed into the Bulaevskii-Ginzburg (narrow-wall) regime.

Protocol:
  1. Build cuboid, define couplings, set anisotropy.
  2. Pin the end slabs and seed a two-domain (+z | -z) initial configuration.
  3. Burn in at T_hot, anneal T_hot -> T_sample (pure Metropolis, static field).
  4. Sample the x-resolved magnetisation profile <S>(x) and net magnetisation.

Output (HDF5):
  /domain_wall/x        [nx]     integer x-slab coordinate (nx = 2*Lx)
  /domain_wall/M        [nx,3]   slab-averaged <S> profile (column 2 = Mz = OP)
  /domain_wall/M_net    [3]      net magnetisation per spin
  /domain_wall/D, T_sample, Lx, Ly, Lz, pin_thickness, B[3], easy_axis[3]
  /energy                        annealing energy trace
  /geometry                      lattice metadata
*/

using namespace std;
using namespace CMC;


int main (int argc, char *argv[]) {
    ///////////////////////////////////////////////////////////////////////////
    /// Setup and CLI
    argparse::ArgumentParser prog("DC");

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
        .help("Temperature to begin annealing from (default: 40)")
        .default_value(40.0)
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
        .help("Number of annealing steps (temperature grid points)")
        .default_value(static_cast<size_t>(200))
        .scan<'i', size_t>();
    prog.add_argument("--n_burn_in")
        .help("Number of sweeps at T_hot before annealing")
        .default_value(static_cast<size_t>(200))
        .scan<'i', size_t>();
    prog.add_argument("--n_sweep", "-w")
        .help("Number of MC sweeps per anneal step / between profile samples")
        .default_value(static_cast<size_t>(16))
        .scan<'i', size_t>();
    prog.add_argument("--n_sample")
        .help("Number of profile measurements to average at T_sample")
        .default_value(static_cast<size_t>(200))
        .scan<'i', size_t>();

    /// PHYSICS
    declare_LJ123(prog);          // positional L, --J1..--Jc1, -B
    declare_cuboid_and_aniso(prog); // --Lx/--Ly/--Lz, --D_sia, --easy_axis, --pin_thickness

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

    const size_t seed       = prog.get<size_t>("--seed");
    const double T_hot      = prog.get<double>("--T_hot");
    const double T_sample   = prog.get<double>("--T_sample");
    const size_t n_anneal   = prog.get<size_t>("--n_anneal_steps");
    const size_t n_burn_in  = prog.get<size_t>("--n_burn_in");
    const size_t n_sweep    = prog.get<size_t>("--n_sweep");
    const size_t n_sample   = prog.get<size_t>("--n_sample");

    const auto dims = get_cuboid_dims(prog);
    const int Lx = dims[0], Ly = dims[1], Lz = dims[2];
    int pin_thickness = prog.get<int>("--pin_thickness");
    if (pin_thickness < 2)
        fprintf(stderr, "[warn] pin_thickness=%d < 2: bulk boundary may be "
                        "corrupted by wrap-around bonds.\n", pin_thickness);

    const double D = prog.get<double>("--D_sia");
    vector3::vec3d easy_axis = from_std_vector(prog.get<std::vector<double>>("--easy_axis"));

    auto lat = build_cuboid_lat(prog);
    auto mc  = build_mc_run(prog, lat, seed);
    mc.set_anisotropy(D, easy_axis);

    // Number of integer x-slabs, derived from the actual lattice (the supercell
    // matrix is specified in units of unit cells, each 2 integer-x wide, so the
    // integer-x extent is not simply 2*Lx).
    int nx = 0;
    for (const auto& s : lat.get_objects<HeisenbergSpin>())
        nx = std::max(nx, static_cast<int>(s.ipos[0]) + 1);
    if (2 * pin_thickness >= nx)
        throw runtime_error("pin_thickness too large: no bulk spins remain");

    const auto B_vec = prog.get<std::vector<double>>("-B");
    vector3::vec3d B = from_std_vector(B_vec);
    mc.set_global_field(B);

    // Normalised easy axis for building the pinned directions.
    const double axn = std::sqrt(easy_axis[0]*easy_axis[0]
                               + easy_axis[1]*easy_axis[1]
                               + easy_axis[2]*easy_axis[2]);
    vector3::vec3d ehat = (axn > 1e-12)
        ? vector3::vec3d{easy_axis[0]/axn, easy_axis[1]/axn, easy_axis[2]/axn}
        : vector3::vec3d{0,0,1};

    std::stringstream name;
    name << prog.get<std::string>("--prefix") << DELIM
         << name_LJ123(prog)
         << "Lx=" << Lx << ",Ly=" << Ly << ",Lz=" << Lz << DELIM
         << "D="  << D << DELIM
         << "B="  << B_vec[0] << "," << B_vec[1] << "," << B_vec[2] << DELIM
         << "seed=" << seed << DELIM
         << "Tc=" << T_sample << DELIM
         << "sw=" << n_sweep << DELIM;

    ///////////////////////////////////////////////////////////////////////////
    /// Pin end slabs and seed a two-domain configuration

    size_t n_pinned = 0;
    for (auto& s : lat.get_objects<HeisenbergSpin>()) {
        const int x = static_cast<int>(s.ipos[0]);
        if (x < pin_thickness) {
            s.S = ehat;  s.pinned = true;  n_pinned++;
        } else if (x >= nx - pin_thickness) {
            s.S = -1.0 * ehat;  s.pinned = true;  n_pinned++;
        } else {
            // Bulk: sharp step so a single wall seeds near mid-x.
            s.S = (x < nx / 2) ? ehat : (-1.0 * ehat);
        }
    }

    const size_t n_spins = lat.get_objects<HeisenbergSpin>().size();
    printf("DC: Lx=%d Ly=%d Lz=%d  n_spins=%zu  pinned=%zu (thickness=%d)\n",
           Lx, Ly, Lz, n_spins, n_pinned, pin_thickness);
    printf("    D=%.4g  easy_axis=(%.3g,%.3g,%.3g)  B=(%.3g,%.3g,%.3g)\n",
           D, ehat[0], ehat[1], ehat[2], B[0], B[1], B[2]);
    printf("    T_hot=%.3g  T_sample=%.3g  n_anneal=%zu  n_burn_in=%zu\n",
           T_hot, T_sample, n_anneal, n_burn_in);

    ///////////////////////////////////////////////////////////////////////////
    /// Output file

    auto file_path = outdir / (name.str() + ".out.h5");
    hid_t file_id = H5Fcreate(file_path.string().c_str(), H5F_ACC_TRUNC,
                               H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0)
        throw std::runtime_error("Failed to create HDF5 file: " + file_path.string());

    energy_manager e_manager;
    auto [T_grid, _] = generate_T_profile(T_hot, T_sample, {}, n_anneal);

    ///////////////////////////////////////////////////////////////////////////
    /// Stage 0: burn-in at T_hot (pure Metropolis — anisotropy active)

    printf("\nBurning in (%zu sweeps at T=%.3g)...\n", n_burn_in, T_hot);
    for (size_t i = 0; i < n_burn_in; i++)
        mc.sweep_metropolis(T_hot);

    ///////////////////////////////////////////////////////////////////////////
    /// Stage 1: anneal T_hot -> T_sample

    printf("Annealing (%zu steps)...\n", T_grid.size());
    for (size_t i = 0; i < T_grid.size(); ++i) {
        const double T = T_grid[i];
        size_t accepted = 0;
        e_manager.new_T(T);
        for (size_t n = 0; n < n_sweep; n++)
            accepted += mc.sweep_metropolis(T);
        e_manager.sample(mc.total_energy_per_unit_cell());
        if (i % 20 == 0 || i + 1 == T_grid.size())
            printf("[anneal] %4zu  T=%.3e  E=%.3e  acc=%.1f%%\n",
                   i, T, e_manager.curr_E(),
                   accepted * 100.0 / n_spins / n_sweep);
    }

    ///////////////////////////////////////////////////////////////////////////
    /// Stage 2: sample the x-resolved magnetisation profile

    printf("\nSampling domain-wall profile (%zu measurements)...\n", n_sample);

    std::vector<double> prof(3 * nx, 0.0);     // summed <S> per slab
    std::vector<double> prof_cnt(nx, 0.0);      // accumulated spin count per slab
    double m_net[3] = {0, 0, 0};

    for (size_t m = 0; m < n_sample; m++) {
        for (size_t n = 0; n < n_sweep; n++)
            mc.sweep_metropolis(T_sample);

        for (const auto& s : lat.get_objects<HeisenbergSpin>()) {
            const int x = static_cast<int>(s.ipos[0]);
            prof[3*x + 0] += s.S[0];
            prof[3*x + 1] += s.S[1];
            prof[3*x + 2] += s.S[2];
            prof_cnt[x]   += 1.0;
            m_net[0] += s.S[0];
            m_net[1] += s.S[1];
            m_net[2] += s.S[2];
        }
    }

    // Normalise: per-slab average and net magnetisation per spin.
    for (int x = 0; x < nx; x++) {
        const double c = (prof_cnt[x] > 0) ? prof_cnt[x] : 1.0;
        prof[3*x+0] /= c;  prof[3*x+1] /= c;  prof[3*x+2] /= c;
    }
    const double inv = 1.0 / (static_cast<double>(n_spins) * static_cast<double>(n_sample));
    m_net[0] *= inv;  m_net[1] *= inv;  m_net[2] *= inv;

    printf("M_net = (%.4e, %.4e, %.4e)\n", m_net[0], m_net[1], m_net[2]);

    ///////////////////////////////////////////////////////////////////////////
    /// Write /domain_wall

    {
        hid_t grp = H5Gcreate2(file_id, "/domain_wall",
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (grp < 0) throw std::runtime_error("Failed to create /domain_wall group");

        auto write_scalar_d = [&](const char* n_, double v) {
            hid_t sp = H5Screate(H5S_SCALAR);
            hid_t ds = H5Dcreate2(grp, n_, H5T_NATIVE_DOUBLE, sp,
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v);
            H5Dclose(ds); H5Sclose(sp);
        };
        auto write_scalar_i = [&](const char* n_, int v) {
            hid_t sp = H5Screate(H5S_SCALAR);
            hid_t ds = H5Dcreate2(grp, n_, H5T_NATIVE_INT, sp,
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dwrite(ds, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v);
            H5Dclose(ds); H5Sclose(sp);
        };
        auto write_vec3 = [&](const char* n_, const double v[3]) {
            const hsize_t d = 3;
            hid_t sp = H5Screate_simple(1, &d, nullptr);
            hid_t ds = H5Dcreate2(grp, n_, H5T_NATIVE_DOUBLE, sp,
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v);
            H5Dclose(ds); H5Sclose(sp);
        };

        // x-slab coordinates
        std::vector<int32_t> xs(nx);
        for (int x = 0; x < nx; x++) xs[x] = x;
        {
            const hsize_t d = static_cast<hsize_t>(nx);
            hid_t sp = H5Screate_simple(1, &d, nullptr);
            hid_t ds = H5Dcreate2(grp, "x", H5T_STD_I32LE, sp,
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dwrite(ds, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, xs.data());
            H5Dclose(ds); H5Sclose(sp);
        }
        // Profile M[nx,3]
        {
            const hsize_t d2[2] = {static_cast<hsize_t>(nx), 3};
            hid_t sp = H5Screate_simple(2, d2, nullptr);
            hid_t ds = H5Dcreate2(grp, "M", H5T_NATIVE_DOUBLE, sp,
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, prof.data());
            H5Dclose(ds); H5Sclose(sp);
        }

        const double ax_arr[3] = {ehat[0], ehat[1], ehat[2]};
        const double B_arr[3]  = {B[0], B[1], B[2]};
        write_vec3("M_net", m_net);
        write_vec3("easy_axis", ax_arr);
        write_vec3("B", B_arr);
        write_scalar_d("D", D);
        write_scalar_d("T_sample", T_sample);
        write_scalar_i("Lx", Lx);
        write_scalar_i("Ly", Ly);
        write_scalar_i("Lz", Lz);
        write_scalar_i("pin_thickness", pin_thickness);

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
