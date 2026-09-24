#include "KarstNSim/compact_sampler_grid.h"

#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

int main() {
    try {
        constexpr int dx = 43, dy = 25, dz = 17;
        KarstNSim::CompactSamplerGrid grid(dx, dy, dz);
        struct Point { int x, y, z, r; };
        std::vector<Point> points;
        std::mt19937 rng(138);
        for (int i = 0; i < 300; ++i) {
            Point p{static_cast<int>(rng() % dx), static_cast<int>(rng() % dy),
                static_cast<int>(rng() % dz), static_cast<int>(rng() % 21)};
            points.push_back(p);
            grid.insert(p.x, p.y, p.z, p.r);
        }
        // Every cell gets the same membership as the original replicated-list
        // sampler, including clipped cubes, radius zero and overlapping points.
        for (int z = 0; z < dz; ++z) {
            for (int y = 0; y < dy; ++y) {
                for (int x = 0; x < dx; ++x) {
                    std::vector<std::int32_t> expected, actual;
                    for (std::size_t i = 0; i < points.size(); ++i) {
                        const auto p = points[i];
                        if (std::abs(p.x - x) <= p.r && std::abs(p.y - y) <= p.r && std::abs(p.z - z) <= p.r)
                            expected.push_back(static_cast<std::int32_t>(i));
                    }
                    grid.any_covering(x, y, z, grid.max_rho(), [&](std::int32_t id) {
                        actual.push_back(id);
                        return false;
                    });
                    std::sort(actual.begin(), actual.end());
                    if (expected != actual) throw std::runtime_error("sampler membership changed");
                    const bool hit = grid.any_covering(x, y, z, grid.max_rho(), [](std::int32_t) { return true; });
                    if (hit != !expected.empty()) throw std::runtime_error("sampler early-out changed");
                }
            }
        }
        // Query arithmetic remains defined even for a radius spanning the domain.
        if (!grid.any_covering(0, 0, 0, std::numeric_limits<int>::max(), [](std::int32_t) { return true; }))
            throw std::runtime_error("large query radius missed a covering point");
        bool threw = false;
        try { KarstNSim::CompactSamplerGrid invalid(0, 1, 1); }
        catch (const std::invalid_argument&) { threw = true; }
        if (!threw) throw std::runtime_error("invalid dimensions accepted");
        std::cout << "Compact sampler matches original cell memberships.\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
