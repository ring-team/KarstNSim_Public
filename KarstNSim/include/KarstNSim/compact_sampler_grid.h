#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <vector>

namespace KarstNSim {

// The original sampler inserts each point into every cell in its scan_area2
// cube. This index stores the point once in a coarse bucket, then filters by
// that same cube during lookup. Iteration order is immaterial: the caller only
// tests whether any covering point is too close, without changing RNG state.
class CompactSamplerGrid {
public:
    // Four original cells per bucket edge balances storage and query work.
    explicit CompactSamplerGrid(int dx, int dy, int dz) : dx_(dx), dy_(dy), dz_(dz) {
        if (dx <= 0 || dy <= 0 || dz <= 0)
            throw std::invalid_argument("Poisson grid dimensions must be positive.");
        bx_ = (dx - 1) / bucket_width + 1;
        by_ = (dy - 1) / bucket_width + 1;
        const std::size_t bz = (dz - 1) / bucket_width + 1;
        const auto limit = heads_.max_size();
        if (static_cast<std::size_t>(bx_) > limit / static_cast<std::size_t>(by_) ||
            static_cast<std::size_t>(bx_) * by_ > limit / bz)
            throw std::length_error("Poisson bucket grid is too large.");
        heads_.assign(static_cast<std::size_t>(bx_) * by_ * bz, -1);
    }

    void insert(int x, int y, int z, int radius) {
        check_cell(x, y, z);
        if (radius < 0) throw std::invalid_argument("Poisson cell radius must not be negative.");
        if (points_.size() >= static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()))
            throw std::length_error("Too many Poisson sampling points.");
        const std::size_t b = bucket(x / bucket_width, y / bucket_width, z / bucket_width);
        const auto id = static_cast<std::int32_t>(points_.size());
        points_.push_back({x, y, z, radius, heads_[b]});
        heads_[b] = id;
        max_radius_ = std::max(max_radius_, radius);
    }

    int max_rho() const { return max_radius_; }

    // q=max_rho() returns the full original membership set. A smaller q can
    // be used only when the caller proves that omitted points cannot reject
    // its candidate. IDs are insertion order, zero based.
    template<class Predicate>
    bool any_covering(int x, int y, int z, int q, Predicate&& predicate) const {
        check_cell(x, y, z);
        if (q < 0) throw std::invalid_argument("Poisson query radius must not be negative.");
        const int x0 = lower(x, q), x1 = upper(x, q, dx_);
        const int y0 = lower(y, q), y1 = upper(y, q, dy_);
        const int z0 = lower(z, q), z1 = upper(z, q, dz_);
        for (int bz = z0; bz <= z1; ++bz)
            for (int by = y0; by <= y1; ++by)
                for (int bx = x0; bx <= x1; ++bx)
                    for (std::int32_t id = heads_[bucket(bx, by, bz)]; id >= 0; id = points_[id].next) {
                        const Entry& p = points_[id];
                        // Coordinates are nonnegative int32, so subtraction is safe.
                        if (std::abs(p.x - x) <= p.radius && std::abs(p.y - y) <= p.radius &&
                            std::abs(p.z - z) <= p.radius && predicate(id))
                            return true;
                    }
        return false;
    }

    std::size_t storage_bytes() const {
        return heads_.capacity() * sizeof(std::int32_t) + points_.capacity() * sizeof(Entry);
    }

private:
    static constexpr int bucket_width = 4;
    struct Entry { std::int32_t x, y, z, radius, next; };
    static int lower(int c, int q) { return std::max<std::int64_t>(static_cast<std::int64_t>(c) - q, 0) / bucket_width; }
    static int upper(int c, int q, int dim) {
        return std::min<std::int64_t>(static_cast<std::int64_t>(c) + q, dim - 1) / bucket_width;
    }
    void check_cell(int x, int y, int z) const {
        if (x < 0 || x >= dx_ || y < 0 || y >= dy_ || z < 0 || z >= dz_)
            throw std::out_of_range("Poisson cell is outside the sampling grid.");
    }
    std::size_t bucket(int x, int y, int z) const {
        return static_cast<std::size_t>(x) + static_cast<std::size_t>(bx_) *
            (static_cast<std::size_t>(y) + static_cast<std::size_t>(by_) * z);
    }
    int dx_, dy_, dz_, bx_ = 0, by_ = 0, max_radius_ = 0;
    std::vector<std::int32_t> heads_;
    std::vector<Entry> points_;
};

} // namespace KarstNSim
