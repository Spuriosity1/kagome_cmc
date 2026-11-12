#include <complex>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include "H5Apublic.h"
#include "H5Gpublic.h"
#include "H5Ppublic.h"
#include "H5Tpublic.h"

#include "cli_bits.hpp"
#include "geometry/kagome.hpp"
#include "format_bits.hpp"

// LatticeIndexing lives in supercell.hpp (included via cli_bits → kagome → MC.hpp).
// We use it directly; no Supercell<HeisenbergSpin> is ever constructed.

using namespace std;

// Luttinger–Tisza analysis for the kagome spin model.
//
// The exchange couplings are in general anisotropic 3×3 matrices (the NN/J1
// bond carries a Dzyaloshinskii–Moriya part; the rest are isotropic Heisenberg),
// so the Fourier-transformed coupling matrix M(k) is (3·n_sl)×(3·n_sl) complex
// Hermitian: block M[3α:3α+3, 3β:3β+3] accumulates J_bond·e^{ik·d} over directed
// bonds α→β. The LT lower bound on the energy per spin is ½·min_k λ_min(M(k)).
//
// Each undirected bond is listed once in kagome::vecs (from one end only), so
// build_bonds emits it twice: the directed bond α→β and its Hermitian partner
// β→α (displacement −d, matrix Jᵀ). That makes M(k) Hermitian by construction
// and counts every undirected bond twice, which is what the ½ gives back.
//
// Bond orientation (DM sign): the directed bond from the higher-sublattice spin
// to the lower one carries the matrix M, the reverse carries Mᵀ. This matches
// the MC convention H_bond = S_hiᵀ M S_lo (set by the sublattice-major memory
// layout; see cli_bits::build_mc_run), so ltgs and the MC share a Hamiltonian.

// Precomputed directed bond for the LT matrix.
struct LTBond {
    int             alpha, beta;  // source / target sublattice (0 .. n_sl-1)
    ipos_t          disp;         // full integer displacement from alpha to beta
    Eigen::Matrix3d J;            // 3×3 coupling matrix for this directed bond
};

static vector<LTBond> build_bonds(
        LatticeIndexing& lat,
        const vector<ipos_t>& sl_pos,
        const std::map<std::string, double>& J,
        double D_dm, double d_dm)
{
    const int n_sl = static_cast<int>(sl_pos.size());

    // Canonical (higher-sl → lower-sl) 3×3 matrix for each active coupling.
    std::vector<std::pair<std::string, Eigen::Matrix3d>> specs;

    const double J1val = J.count("J1") ? J.at("J1") : 0.0;
    const double dD    = d_dm * D_dm;
    const bool dm_active = (std::abs(D_dm) > 1e-12 || std::abs(dD) > 1e-12);
    if (dm_active || std::abs(J1val) > 1e-10) {
        Eigen::Matrix3d M;
        M <<  J1val,  dD,     D_dm,
             -dD,     J1val,  0.0,
             -D_dm,   0.0,    J1val;
        specs.emplace_back("J1", M);
    }
    for (const auto& [cn, v] : J) {
        if (cn == "J1") continue;
        if (std::abs(v) > 1e-10)
            specs.emplace_back(cn, v * Eigen::Matrix3d::Identity());
    }

    vector<LTBond> bonds;
    for (const auto& [cn, Mc] : specs) {
        const auto& dist = kagome::vecs.at(cn); // per-sublattice displacement table
        for (int sl = 0; sl < n_sl; sl++) {
            for (const auto& v : dist[sl]) {
                ipos_t ref = sl_pos[sl] + v;
                lat.get_supercell_IDX(ref); // wrap ref into the primitive cell
                int beta = -1;
                for (int s = 0; s < n_sl; s++) {
                    if (sl_pos[s] == ref) { beta = s; break; }
                }
                assert(beta >= 0 && "bond target not found in sl_positions");
                // Higher-sl source → carry M; lower-sl source → carry Mᵀ.
                Eigen::Matrix3d Jb = (sl > beta) ? Mc : Mc.transpose();
                // Emit the directed bond α→β together with its Hermitian partner
                // β→α (displacement −v, matrix Jbᵀ). This makes M(k) Hermitian by
                // construction (block[β,α] = block[α,β]†) and counts every
                // undirected bond twice, which is what the ½ in E_per_spin assumes.
                ipos_t vneg{-v[0], -v[1], -v[2]};
                bonds.push_back({sl, beta, v, Jb});
                bonds.push_back({beta, sl, vneg, Jb.transpose()});
            }
        }
    }
    return bonds;
}


int main(int argc, char* argv[])
{
    argparse::ArgumentParser prog("ltgs");

    prog.add_argument("--output_dir", "-o")
        .help("Path to output directory");

    declare_LJ123(prog);
    declare_DM(prog);

    try {
        prog.parse_args(argc, argv);
    } catch (const std::exception& err) {
        cerr << err.what() << "\n" << prog;
        return 1;
    }

    filesystem::path outdir;
    if (prog.is_used("--output_dir")){
        const string outdir_s = prog.get<string>("output_dir");
        outdir = outdir_s;
        if (!filesystem::exists(outdir))
            throw runtime_error("output_dir does not exist: " + outdir_s);
    }

    const std::map<std::string, double> J = get_J(prog);
    const double D_dm = prog.get<double>("--D");
    const double d_dm = prog.get<double>("--d");

    // Build LatticeIndexing directly — no Supercell<HeisenbergSpin> needed.
    // Primitive cell matches KagomeLayer(); a diagonal L×L×L supercell (in
    // primitive-cell units) gives a clean L³ k-grid for the BZ sweep.
    const int L = prog.get<int>("L");
    const imat33_t prim_cell = imat33_t::from_cols({2,0,0},{0,2,0},{0,0,1});
    const imat33_t supercell = imat33_t::from_cols({L,0,0},{0,L,0},{0,0,L});
    LatticeIndexing lat(prim_cell, supercell);

    // Sublattice positions (kagome: 3 sites per cell). Wrap each into the
    // primitive cell so positions match what get_supercell_IDX returns. The
    // index order (0,1,2) matches KagomeLayer()'s pyro_sl, hence the MC's
    // sublattice-major layout — so "higher-sl" means the same thing here.
    vector<ipos_t> sl_pos;
    for (const auto& p : std::vector<ipos_t>{ {0,0,0}, {1,0,0}, {1,1,0} }) {
        ipos_t pos = p;
        lat.wrap_primitive(pos);
        sl_pos.push_back(pos);
    }
    const int n_sl = static_cast<int>(sl_pos.size()); // = 3

    const int Nk           = lat.num_primitive_cells();
    const ivec3_t k_dims   = lat.size();
    const auto bonds       = build_bonds(lat, sl_pos, J, D_dm, d_dm);

    // -----------------------------------------------------------------------
    // Sweep BZ: build M(k) and track the minimum eigenvalue.
    //
    // M(k) is (3·n_sl)×(3·n_sl) complex Hermitian; block M[3α:3α+3, 3β:3β+3]
    // accumulates J_bond·e^{ik·d}. Hermiticity follows from the reverse bond
    // (β→α, disp=−d, matrix Jᵀ) contributing the conjugate block.
    // -----------------------------------------------------------------------
    using MatC = Eigen::MatrixXcd;
    const int M_dim = 3 * n_sl;

    vector<double> eigenvalue_map(Nk);
    double E_min = numeric_limits<double>::max();
    idx3_t k_star_idx{};
    Eigen::VectorXcd eigvec_star(M_dim);

    for (int k_flat = 0; k_flat < Nk; k_flat++) {
        const idx3_t Q   = lat.idx3_from_flat(k_flat);
        const auto k_vec = lat.wavevector_from_idx3(Q);

        MatC M = MatC::Zero(M_dim, M_dim);
        for (const auto& b : bonds) {
            const double phase =
                k_vec[0] * (double)b.disp[0] +
                k_vec[1] * (double)b.disp[1] +
                k_vec[2] * (double)b.disp[2];
            complex<double> eph(cos(phase), sin(phase));
            M.block<3,3>(3*b.alpha, 3*b.beta) +=
                b.J.cast<complex<double>>() * eph;
        }
        // M is Hermitian by construction (full directed-bond sweep)

        Eigen::SelfAdjointEigenSolver<MatC> eigs(M);
        const double lmin = eigs.eigenvalues()(0);
        eigenvalue_map[k_flat] = lmin;

        if (lmin < E_min) {
            E_min       = lmin;
            k_star_idx  = Q;
            eigvec_star = eigs.eigenvectors().col(0);
        }
    }

    // Energy per spin: each undirected bond is counted twice in M(k) → ½.
    const double E_per_spin = 0.5 * E_min;

    const auto k_star_vec = lat.wavevector_from_idx3(k_star_idx);

    printf("LT minimum: λ_min = %.6f  E/spin = %.6f\n", E_min, E_per_spin);
    // Per-sublattice weight: sum of |components|² over the 3 spin dofs.
    printf("Eigvec sublattice |s_α|² :");
    for (int s = 0; s < n_sl; s++) {
        double w = norm(eigvec_star(3*s)) + norm(eigvec_star(3*s+1)) + norm(eigvec_star(3*s+2));
        printf(" %.4f", w);
    }
    printf("\n");

    printf("k*  idx = [%lld, %lld, %lld]\n",
           (long long)k_star_idx[0],
           (long long)k_star_idx[1],
           (long long)k_star_idx[2]);
    printf("k*  vec = [%.4f, %.4f, %.4f]  (rad / coord-unit)\n",
           k_star_vec[0], k_star_vec[1], k_star_vec[2]);


    // -----------------------------------------------------------------------
    // HDF5 output
    // -----------------------------------------------------------------------

    // exit if no output file provided
    if (! prog.is_used("--output_dir")) return 0;

    stringstream name;
    name << name_LJ123(prog) << name_DM(prog);
    auto file_path = outdir / (name.str() + "lt.out.h5");

    hid_t fid = H5Fcreate(file_path.string().c_str(),
                           H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (fid < 0)
        throw runtime_error("Failed to create " + file_path.string());

    auto write_1d = [&](const char* dname, hid_t type,
                        hsize_t len, const void* data) {
        hid_t sp = H5Screate_simple(1, &len, nullptr);
        hid_t ds = H5Dcreate2(fid, dname, type, sp,
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
        H5Dclose(ds);
        H5Sclose(sp);
    };

    // eigenvalue_map [Nk]
    write_1d("eigenvalue_map", H5T_NATIVE_DOUBLE,
             static_cast<hsize_t>(Nk), eigenvalue_map.data());

    // k_star_idx [3]
    {
        int64_t ksi[3] = { k_star_idx[0], k_star_idx[1], k_star_idx[2] };
        write_1d("k_star_idx", H5T_NATIVE_INT64, 3, ksi);
    }

    // k_star_vec [3]
    {
        double ksv[3] = { k_star_vec[0], k_star_vec[1], k_star_vec[2] };
        write_1d("k_star_vec", H5T_NATIVE_DOUBLE, 3, ksv);
    }

    // E_min scalar (energy per spin)
    {
        hid_t sp = H5Screate(H5S_SCALAR);
        hid_t ds = H5Dcreate2(fid, "E_min", H5T_NATIVE_DOUBLE, sp,
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &E_per_spin);
        H5Dclose(ds);
        H5Sclose(sp);
    }

    // eigenvector [3*n_sl, 2] (re/im; rows ordered [spin-x,y,z for sl-0, sl-1, ...])
    {
        vector<double> ev(2 * M_dim);
        for (int s = 0; s < M_dim; s++) {
            ev[2*s]   = eigvec_star(s).real();
            ev[2*s+1] = eigvec_star(s).imag();
        }
        hsize_t dims[2] = { static_cast<hsize_t>(M_dim), 2 };
        hid_t sp = H5Screate_simple(2, dims, nullptr);
        hid_t ds = H5Dcreate2(fid, "eigenvector", H5T_NATIVE_DOUBLE, sp,
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ev.data());
        H5Dclose(ds);
        H5Sclose(sp);
    }

    // sl_positions [n_sl, 3]
    {
        vector<int64_t> sp_buf(3 * n_sl);
        for (int s = 0; s < n_sl; s++)
            for (int a = 0; a < 3; a++) sp_buf[3*s + a] = sl_pos[s][a];
        hsize_t dims[2] = { static_cast<hsize_t>(n_sl), 3 };
        hid_t sp = H5Screate_simple(2, dims, nullptr);
        hid_t ds = H5Dcreate2(fid, "sl_positions", H5T_NATIVE_INT64, sp,
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(ds, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, sp_buf.data());
        H5Dclose(ds);
        H5Sclose(sp);
    }

    // /geometry group: lattice and reciprocal vectors
    {
        hid_t grp = H5Gcreate2(fid, "geometry",
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        auto write_mat = [&](const char* dname, hid_t type, const auto& matrix) {
            hsize_t dims[2] = {3, 3};
            decltype(matrix(0,0)) buf[9];
            for (int i = 0; i < 9; i++) buf[i] = matrix(i/3, i%3);
            hid_t sp = H5Screate_simple(2, dims, nullptr);
            hid_t ds = H5Dcreate2(grp, dname, type, sp,
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dwrite(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
            H5Dclose(ds);
            H5Sclose(sp);
        };
        write_mat("index_cell",    H5T_NATIVE_INT64,  lat.get_lattice_vectors());
        write_mat("recip_vectors", H5T_NATIVE_DOUBLE, lat.get_reciprocal_lattice_vectors());

        H5Gclose(grp);
    }

    // k_dims attribute on root
    {
        hsize_t adim = 3;
        int kd[3] = { (int)k_dims[0], (int)k_dims[1], (int)k_dims[2] };
        hid_t sp = H5Screate_simple(1, &adim, nullptr);
        hid_t at = H5Acreate2(fid, "k_dims", H5T_NATIVE_INT, sp,
                               H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(at, H5T_NATIVE_INT, kd);
        H5Aclose(at);
        H5Sclose(sp);
    }

    H5Fclose(fid);
    cout << "Saved to " << file_path << "\n";

    return 0;
}
