#pragma once

#include "MC.hpp"
#include "abstract_manager.hpp"

#include "H5Gpublic.h"
#include "H5Ipublic.h"
#include "H5Ppublic.h"
#include "H5Spublic.h"
#include "H5Dpublic.h"
#include "H5Tpublic.h"

#include <array>
#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <vector>


// Per-temperature order-parameter accumulator.
//
// Measures, each time sample() is called (wired to fire whenever the energy is
// measured in anneal.cpp):
//
//   * the scalar spin chirality, averaged over all registered triangles, and
//   * the sublattice-resolved magnetisation m_mu^alpha = <S^alpha>_{cells},
//     i.e. the per-sublattice average of each spin component over all cells.
//
// As with energy_manager, raw sums (not means) are accumulated on disk together
// with n_samples, so postprocessing divides by n_samples to get expectation
// values and by (sum of squares) to get the fluctuations.
//
// HDF5 output layout (write_group, default group "/order"):
//   scalar_chirality    [n_T]            sum over samples of <chi>_triangles
//   scalar_chirality_2  [n_T]            sum over samples of <chi>^2
//   magnetization       [n_T, n_sl, 3]   sum over samples of m_mu^alpha
//   magnetization_2     [n_T, n_sl, 3]   sum over samples of (m_mu^alpha)^2
//   T_list              [n_T]
//   n_samples           [n_T]
//   sl_positions        [n_sl, 3]        int64 sublattice reference-cell coords
class order_manager : public abstract_manager {
    using Spin = CMC::HeisenbergSpin;

    std::vector<std::array<const vector3::vec3d*, 3>> triangles;

    std::vector<double> scalar_spin_chirality;
    std::vector<double> scalar_spin_chirality2;

    // Flat per-temperature magnetisation accumulators, stride n_sl_*3 per
    // temperature. Index: ((t * n_sl_) + sl) * 3 + component.
    std::vector<double> magnetization;
    std::vector<double> magnetization2;

    CMC::Lattice* sc;
    int n_sl_;   // number of sublattices
    int n_prim_; // number of primitive cells (spins per sublattice)

    void on_new_temp() override {
        scalar_spin_chirality.push_back(0);
        scalar_spin_chirality2.push_back(0);
        magnetization.insert(magnetization.end(),
                static_cast<size_t>(n_sl_) * 3, 0.0);
        magnetization2.insert(magnetization2.end(),
                static_cast<size_t>(n_sl_) * 3, 0.0);
    }

    void register_triangles() {
        triangles.clear();
        for (const auto& [I, cell] : sc->enumerate_cells()) {
            std::array<const vector3::vec3d*, 3> tmp;
            for (const auto& [sl, S] : cell.enumerate_objects<Spin>()) {
                tmp[sl] = &(S->S);
            }
            triangles.push_back(std::move(tmp));
        }
    }

    void cache_dims() {
        n_sl_   = static_cast<int>(sc->num_sl<Spin>());
        n_prim_ = sc->lattice.num_primitive_cells();
    }

    public:
    order_manager(CMC::Lattice& sc_) : abstract_manager(), sc(&sc_) {
        cache_dims();
        // watch all the relevant triangles
        register_triangles();
    }

    void point_at(CMC::Lattice& sc2_) {
        sc = &sc2_;
        cache_dims();
        register_triangles();
    }

    void sample() {
        assert(!T_list.empty());

        // --- scalar spin chirality, averaged over triangles ---
        double ssc = 0;
        for (auto T : triangles) {
            using namespace vector3;
            const auto& s0 = *T[0];
            const auto& s1 = *T[1];
            const auto& s2 = *T[2];

            ssc += s0[0] * (s1[1]*s2[2]-s1[2]*s2[1]) +
                s0[1]*(s1[2]*s2[0]-s1[0]*s2[2]) +
                s0[2]*(s1[0]*s2[1]-s1[1]*s2[0]);
        }
        ssc /= triangles.size();

        scalar_spin_chirality[curr_idx]  += ssc;
        scalar_spin_chirality2[curr_idx] += ssc*ssc;

        // --- sublattice-resolved magnetisation m_mu^alpha ---
        // Spins are laid out contiguously by sublattice:
        //   get_objects<Spin>()[sl * n_prim_ + cell].
        auto& spins = sc->get_objects<Spin>();
        for (int sl = 0; sl < n_sl_; ++sl) {
            double m0 = 0, m1 = 0, m2 = 0;
            const int base_spin = sl * n_prim_;
            for (int c = 0; c < n_prim_; ++c) {
                const auto& S = spins[base_spin + c].S;
                m0 += S[0];
                m1 += S[1];
                m2 += S[2];
            }
            m0 /= n_prim_;
            m1 /= n_prim_;
            m2 /= n_prim_;

            const size_t base = (static_cast<size_t>(curr_idx) * n_sl_ + sl) * 3;
            magnetization[base + 0]  += m0;
            magnetization[base + 1]  += m1;
            magnetization[base + 2]  += m2;
            magnetization2[base + 0] += m0*m0;
            magnetization2[base + 1] += m1*m1;
            magnetization2[base + 2] += m2*m2;
        }

        n_samples[curr_idx]++;
    }

    void save(const std::filesystem::path& file_path);
    void write_group(hid_t file_id, const char* group_name = "/order") override;
};


inline void order_manager::save(const std::filesystem::path& file_path) {
    hid_t file_id = H5Fcreate(file_path.string().c_str(), H5F_ACC_TRUNC,
            H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0) {
        throw std::runtime_error("Failed to create HDF5 file: " + file_path.string());
    }
    write_group(file_id);
    H5Fclose(file_id);
}


inline void order_manager::write_group(hid_t file_id, const char* group_name) {
    hid_t data_group = H5Gcreate2(file_id, group_name,
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (data_group < 0) {
        throw std::runtime_error("Failed to create group");
    }

    auto write_dataset = [&](const char* name, hid_t type,
            int rank, const hsize_t* dims, const void* data) {
        hid_t dataspace_id = H5Screate_simple(rank, dims, nullptr);
        if (dataspace_id < 0) {
            H5Gclose(data_group);
            throw std::runtime_error("Failed to create dataspace");
        }
        hid_t dataset_id = H5Dcreate2(data_group, name, type, dataspace_id,
                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (dataset_id < 0) {
            H5Sclose(dataspace_id);
            H5Gclose(data_group);
            throw std::runtime_error(std::string("Failed to create dataset ") + name);
        }
        herr_t status = H5Dwrite(dataset_id, type, H5S_ALL, H5S_ALL,
                H5P_DEFAULT, data);
        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);
        if (status < 0) {
            H5Gclose(data_group);
            throw std::runtime_error(std::string("Failed to write dataset ") + name);
        }
    };

    const hsize_t n_T = T_list.size();

    // Per-temperature scalars.
    const hsize_t dims_T[1] = { n_T };
    write_dataset("scalar_chirality",   H5T_NATIVE_DOUBLE, 1, dims_T,
            scalar_spin_chirality.data());
    write_dataset("scalar_chirality_2", H5T_NATIVE_DOUBLE, 1, dims_T,
            scalar_spin_chirality2.data());
    write_dataset("T_list",             H5T_NATIVE_DOUBLE, 1, dims_T,
            T_list.data());
    write_dataset("n_samples",          H5T_NATIVE_ULLONG, 1, dims_T,
            n_samples.data());

    // Sublattice-resolved magnetisation [n_T, n_sl, 3].
    const hsize_t dims_m[3] = { n_T, static_cast<hsize_t>(n_sl_), 3 };
    write_dataset("magnetization",   H5T_NATIVE_DOUBLE, 3, dims_m,
            magnetization.data());
    write_dataset("magnetization_2", H5T_NATIVE_DOUBLE, 3, dims_m,
            magnetization2.data());

    // Sublattice reference-cell positions [n_sl, 3].
    {
        const auto& sl_pos = std::get<SlPos<Spin>>(sc->sl_positions);
        std::vector<int64_t> flat;
        flat.reserve(sl_pos.size() * 3);
        for (const auto& p : sl_pos) {
            flat.push_back(p[0]);
            flat.push_back(p[1]);
            flat.push_back(p[2]);
        }
        const hsize_t dims_sl[2] = { static_cast<hsize_t>(n_sl_), 3 };
        write_dataset("sl_positions", H5T_NATIVE_INT64, 2, dims_sl, flat.data());
    }

    H5Gclose(data_group);
}
