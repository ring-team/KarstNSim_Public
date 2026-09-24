#pragma once

// Example factory for in-memory KarstNSim requests, as a region worker would build them:
// a synthetic box domain, a synthetic topographic surface, sinks and springs placed in
// the region, automatic connectivity (no connectivity_matrix.txt) and every filesystem
// export disabled. Everything is built in memory; nothing is read from disk.

#include "KarstNSim/library.h"

#include <cmath>
#include <vector>

namespace KarstNSimTest {

    struct RegionRequest {
        float size_x = 400.0f; //!< Region extent along x (model units).
        float size_y = 400.0f; //!< Region extent along y.
        float size_z = 200.0f; //!< Region thickness.
        float origin_x = 0.0f; //!< Region origin (lets neighbouring regions share one world frame).
        float origin_y = 0.0f;
        float origin_z = 0.0f;
        int seed = 1;
        float poisson_radius = 0.08f; //!< Normalized Poisson radius (fraction of the smallest extent).
        int nghb_count = 16;
        int sink_count = 4;
        int spring_count = 2;
        bool noise = false; //!< Simplex noise on the whole cost graph.
        int cycles = 0; //!< Cycle amplification count (0 disables amplification).
        bool sections = false; //!< Simulate equivalent conduit radii (SGS).
    };

    // Topography: a gently tilted plane covering the region with a margin, triangulated on a grid.
    inline KarstNSim::Surface make_topography(const RegionRequest& r) {
        const int cells = 16;
        const float margin = 0.1f;
        const float x0 = r.origin_x - margin * r.size_x;
        const float y0 = r.origin_y - margin * r.size_y;
        const float dx = r.size_x * (1.0f + 2.0f * margin) / cells;
        const float dy = r.size_y * (1.0f + 2.0f * margin) / cells;
        std::vector<::Vector3> nodes;
        std::vector<KarstNSim::Triangle> triangles;
        for (int j = 0; j <= cells; ++j) {
            for (int i = 0; i <= cells; ++i) {
                const float x = x0 + i * dx;
                const float y = y0 + j * dy;
                const float tilt = 0.04f * r.size_z * ((x - r.origin_x) / r.size_x);
                const float z = r.origin_z + 0.97f * r.size_z - tilt;
                nodes.emplace_back(x, y, z);
            }
        }
        const int stride = cells + 1;
        for (int j = 0; j < cells; ++j) {
            for (int i = 0; i < cells; ++i) {
                const int a = j * stride + i;
                triangles.emplace_back(a, a + 1, a + stride);
                triangles.emplace_back(a + 1, a + stride + 1, a + stride);
            }
        }
        return KarstNSim::Surface(nodes, triangles, "synthetic_topography");
    }

    inline KarstNSim::ParamsSource make_region_params(const RegionRequest& r) {
        KarstNSim::ParamsSource p;
        p.karstic_network_name = "region";
        p.selected_seed = r.seed;
        p.number_of_iterations = 1;

        p.domain = KarstNSim::Box(
            ::Vector3(r.origin_x, r.origin_y, r.origin_z),
            ::Vector3(r.size_x, 0.0f, 0.0f),
            ::Vector3(0.0f, r.size_y, 0.0f),
            ::Vector3(0.0f, 0.0f, r.size_z),
            10, 10, 5);
        p.topo_surface = make_topography(r);

        p.use_sampling_points = false;
        p.use_density_property = false;
        p.poisson_radius = r.poisson_radius;
        p.k_pts = 10;
        p.nghb_count = r.nghb_count;

        // Sinks near the top of the region, springs low on its downstream (+x) side.
        for (int i = 0; i < r.sink_count; ++i) {
            const float fy = (i + 0.5f) / r.sink_count;
            const float fx = 0.15f + 0.25f * static_cast<float>(i % 2);
            p.sinks.emplace_back(r.origin_x + fx * r.size_x, r.origin_y + fy * r.size_y,
                r.origin_z + 0.85f * r.size_z);
            p.propsinksindex.push_back(i + 1);
            p.propsinksorder.push_back(1);
        }
        for (int i = 0; i < r.spring_count; ++i) {
            const float fy = (i + 0.5f) / r.spring_count;
            p.springs.emplace_back(r.origin_x + 0.9f * r.size_x, r.origin_y + fy * r.size_y,
                r.origin_z + 0.25f * r.size_z);
            p.propspringsindex.push_back(i + 1);
            p.propspringssurfindex.push_back(0); // no water table: vadose-only outlet channel
        }
        p.allow_single_outlet_connection = true;
        p.use_user_connectivity_matrix = false; // automatic all-2 connectivity, resolved in memory

        p.fraction_karst_perm = 0.9f;
        p.gamma = 2.0f;

        p.use_amplification = r.cycles > 0;
        p.nb_cycles = r.cycles;
        p.min_distance_amplification = 0.1f * r.size_x;
        p.max_distance_amplification = 0.5f * r.size_x;

        p.use_noise = r.noise;
        p.use_noise_on_all = r.noise;
        p.noise_frequency = 4;
        p.noise_octaves = 2;
        p.noise_weight = 10.0f;

        p.simulate_sections = r.sections;
        p.geostat_params.is_used = r.sections;
        if (r.sections) {
            ::GeostatParams& g = p.geostat_params;
            for (int i = 0; i < 200; ++i) {
                const float t = i / 199.0f;
                g.simulation_distribution.push_back(-2.8f + 3.2f * t * t);
            }
            g.global_vario_range = 50; g.global_range_of_neighborhood = 150;
            g.global_vario_sill = 0.92f; g.global_vario_nugget = 0.33f; g.global_vario_model = "Spherical";
            g.interbranch_vario_range = 30; g.interbranch_range_of_neighborhood = 80;
            g.interbranch_vario_sill = 0.92f; g.interbranch_vario_nugget = 0.34f; g.interbranch_vario_model = "Gaussian";
            g.intrabranch_vario_range = 30; g.intrabranch_range_of_neighborhood = 45;
            g.intrabranch_vario_sill = 0.55f; g.intrabranch_vario_nugget = 0.31f; g.intrabranch_vario_model = "Exponential";
            g.number_max_of_neighborhood_points = 16;
            g.nb_points_interbranch = 7;
            g.proportion_interbranch = 0.1f;
        }

        // Filesystem-only exports must be off for run_simulation_memory.
        p.create_vset_sampling = false;
        p.create_nghb_graph = false;
        p.create_nghb_graph_property = false;
        p.create_grid = false;
        p.create_solved_connectivity_matrix = true; // returned in memory
        return p;
    }
}
