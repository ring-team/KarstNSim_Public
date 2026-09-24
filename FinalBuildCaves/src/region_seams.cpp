// Shared generation settings and canonical region seams.
//
// A seam belongs to the canonical neighbour pair (lower key, axis). Both regions
// derive it from that pair plus the shared settings only, so the key, width,
// height and plane position agree bit-for-bit regardless of which region is
// generated first, on which thread, or whether the other region exists at all.

#include "region_internal.hpp"

#include <algorithm>
#include <limits>

namespace fbs::caves::detail {

Shared shared_from_request(const Request& r) {
    Shared s;
    s.world_seed = r.world_seed;
    s.style_id = r.style_id;
    s.size = r.size;
    s.chamber_count = r.chamber_count;
    s.loop_count = r.loop_count;
    s.neighbors = r.neighbors;
    s.chamber_radius_min = r.chamber_radius_min;
    s.chamber_radius_max = r.chamber_radius_max;
    s.chamber_height = r.chamber_height;
    s.tunnel_width = r.tunnel_width;
    s.tunnel_height = r.tunnel_height;
    s.maximum_slope_degrees = r.maximum_slope_degrees;
    s.sampling_radius = r.sampling_radius;
    s.shelves = r.shelves;
    return s;
}

std::vector<std::pair<std::string, std::string>> shared_parameters(const Shared& s) {
    return {
        {"fbs.algorithm", algorithm_version},
        {"fbs.schema_version", std::to_string(schema_version)},
        {"fbs.world_seed", std::to_string(s.world_seed)},
        {"fbs.style_id", std::to_string(s.style_id)},
        {"fbs.size_x", format_double(s.size.x)},
        {"fbs.size_y", format_double(s.size.y)},
        {"fbs.size_z", format_double(s.size.z)},
        {"fbs.chamber_count", std::to_string(s.chamber_count)},
        {"fbs.loop_count", std::to_string(s.loop_count)},
        {"fbs.neighbors", std::to_string(s.neighbors)},
        {"fbs.chamber_radius_min", format_double(s.chamber_radius_min)},
        {"fbs.chamber_radius_max", format_double(s.chamber_radius_max)},
        {"fbs.chamber_height", format_double(s.chamber_height)},
        {"fbs.tunnel_width", format_double(s.tunnel_width)},
        {"fbs.tunnel_height", format_double(s.tunnel_height)},
        {"fbs.maximum_slope_degrees", format_double(s.maximum_slope_degrees)},
        {"fbs.sampling_radius", format_double(s.sampling_radius)},
        {"fbs.shelves", s.shelves ? "1" : "0"},
    };
}

bool shared_from_parameters(const std::map<std::string, std::string>& p, Shared& s) {
    auto get = [&](const char* k, std::string& out) {
        const auto it = p.find(k);
        if (it == p.end()) return false;
        out = it->second;
        return true;
    };
    std::string v;
    std::uint64_t u = 0;
    auto u32 = [&](const char* k, std::uint32_t& out) {
        if (!get(k, v) || !parse_u64(v, u) || u > std::numeric_limits<std::uint32_t>::max()) return false;
        out = static_cast<std::uint32_t>(u);
        return true;
    };
    auto dbl = [&](const char* k, double& out) { return get(k, v) && parse_double(v, out); };
    if (!get("fbs.algorithm", v) || v != algorithm_version) return false;
    if (!get("fbs.schema_version", v) || v != std::to_string(schema_version)) return false;
    if (!get("fbs.world_seed", v) || !parse_u64(v, s.world_seed)) return false;
    if (!get("fbs.style_id", v) || !parse_u64(v, s.style_id)) return false;
    if (!dbl("fbs.size_x", s.size.x) || !dbl("fbs.size_y", s.size.y) || !dbl("fbs.size_z", s.size.z)) return false;
    if (!u32("fbs.chamber_count", s.chamber_count) || !u32("fbs.loop_count", s.loop_count) ||
        !u32("fbs.neighbors", s.neighbors))
        return false;
    if (!dbl("fbs.chamber_radius_min", s.chamber_radius_min) || !dbl("fbs.chamber_radius_max", s.chamber_radius_max) ||
        !dbl("fbs.chamber_height", s.chamber_height) || !dbl("fbs.tunnel_width", s.tunnel_width) ||
        !dbl("fbs.tunnel_height", s.tunnel_height) || !dbl("fbs.maximum_slope_degrees", s.maximum_slope_degrees) ||
        !dbl("fbs.sampling_radius", s.sampling_radius))
        return false;
    if (!get("fbs.shelves", v) || (v != "0" && v != "1")) return false;
    s.shelves = v == "1";
    return true;
}

std::string settings_id(const Shared& s) {
    std::string canonical = "fbs-caves/settings";
    for (const auto& [k, v] : shared_parameters(s)) canonical += "|" + k + "=" + v;
    return stable_id(canonical);
}

std::string region_identity(const Shared& s, const RegionKey& key) {
    return "fbs-caves/1/region/" + settings_id(s) + "/" + std::to_string(key.x) + "/" + std::to_string(key.y) + "/" +
           std::to_string(key.z);
}

Derived derive(const Shared& s) {
    Derived d;
    d.gap = std::max(1.0, 0.25 * s.tunnel_width);
    d.width_max = 1.3 * s.tunnel_width;
    d.height_max = 1.15 * s.tunnel_height;
    d.radius = section_radius(d.width_max, d.height_max);
    d.spacing = std::max({1.0, 0.5 * s.tunnel_width, s.sampling_radius * std::min({s.size.x, s.size.y, s.size.z})});
    d.gateway_radius = std::max(s.chamber_radius_min, 2.0 * d.radius + d.gap + 2.0);
    // Tall enough that every gateway port (seam or internal) sits inside the
    // support-lattice margin whatever the seam elevation.
    d.gateway_height = std::max({1.5 * d.height_max, d.height_max + 2.0, d.radius + 0.5 * d.height_max + 0.5});
    d.corridor = std::max(2.0 * d.radius + d.gap, 0.3 * d.gateway_radius + d.gap + d.radius);
    return d;
}

Face opposite(Face f) {
    const int v = static_cast<int>(f);
    return static_cast<Face>(v % 2 == 0 ? v + 1 : v - 1);
}

bool neighbor_key(const RegionKey& key, Face face, RegionKey& out) {
    out = key;
    std::int64_t* c = face_axis(face) == 0 ? &out.x : face_axis(face) == 1 ? &out.y : &out.z;
    if (face_positive(face)) {
        if (*c == std::numeric_limits<std::int64_t>::max()) return false;
        ++*c;
    } else {
        if (*c == std::numeric_limits<std::int64_t>::min()) return false;
        --*c;
    }
    return true;
}

namespace {
double axis_of(Vec3 v, int a) { return a == 0 ? v.x : a == 1 ? v.y : v.z; }
void set_axis(Vec3& v, int a, double value) { (a == 0 ? v.x : a == 1 ? v.y : v.z) = value; }
} // namespace

Seam derive_seam(const Shared& s, const Derived& d, const RegionKey& key, Face face) {
    RegionKey lower = key;
    if (!face_positive(face)) {
        if (!neighbor_key(key, face, lower))
            throw Error(ErrorCode::invalid_input, std::string("region has no neighbour across face ") + face_label(face));
    } else {
        RegionKey upper;
        if (!neighbor_key(key, face, upper))
            throw Error(ErrorCode::invalid_input, std::string("region has no neighbour across face ") + face_label(face));
    }
    const int axis = face_axis(face);
    const std::string identity = "fbs-caves/1/seam/" + settings_id(s) + "/" + std::to_string(axis) + "/" +
                                 std::to_string(lower.x) + "/" + std::to_string(lower.y) + "/" + std::to_string(lower.z);
    Mixer128 mixer("fbs.caves.seam.v1");
    mixer.str(identity);
    Rng rng(mixer.finish());

    Seam seam;
    seam.key = stable_id(identity);
    seam.width = s.tunnel_width * (1.0 + 0.25 * rng.uniform());
    seam.height = s.tunnel_height * (1.0 + 0.15 * rng.uniform());
    const bool positive = face_positive(face);

    if (axis < 2) {
        const int tangent = axis == 0 ? 1 : 0;
        const double lo = d.radius + d.gap, hi = s.size.z - d.gateway_height - d.gap;
        if (hi < lo)
            throw Error(ErrorCode::constraint, "region height cannot hold a gateway chamber for a horizontal seam");
        const double t = axis_of(s.size, tangent) * (0.25 + 0.5 * rng.uniform());
        const double floor_z = lo + (hi - lo) * rng.uniform();
        Vec3 p;
        set_axis(p, axis, positive ? axis_of(s.size, axis) : 0.0);
        set_axis(p, tangent, t);
        p.z = floor_z + 0.5 * seam.height;
        seam.position = p;
        seam.port = p;
        set_axis(seam.port, axis, positive ? axis_of(s.size, axis) - d.corridor : d.corridor);
        seam.facing = {};
        set_axis(seam.facing, axis, positive ? 1.0 : -1.0);
        seam.slope_degrees = 0;
        return seam;
    }

    // Vertical seam: a straight ramp through the shared horizontal plane.
    const double alpha = std::min(s.maximum_slope_degrees - 2.0 * slope_margin_degrees, 35.0);
    if (!(alpha > 0.5))
        throw Error(ErrorCode::constraint, "vertical face requested but the slope limit leaves no feasible ramp");
    // Plane to seam-port centre on each side: keeps the lower gateway ceiling
    // below the plane and every gateway port inside the lattice margin.
    const double delta = std::max(d.gateway_height - 0.5 * seam.height + d.gap, d.radius + d.gap + 0.5 * seam.height);
    const double run = delta / std::tan(alpha * pi / 180.0);
    const double reach = run + 2.0 * d.gateway_radius + d.gap + d.radius;
    const double upper_floor = delta - 0.5 * seam.height;
    if (2.0 * reach >= std::min(s.size.x, s.size.y) || upper_floor + d.gateway_height + d.gap > s.size.z ||
        s.size.z - delta - 0.5 * seam.height < d.gap)
        throw Error(ErrorCode::constraint,
                    "vertical face infeasible: the slope limit needs a " + format_double(run) +
                        " m horizontal ramp that does not fit the region");
    const double angle = 2.0 * pi * rng.uniform();
    const Vec3 dir{std::cos(angle), std::sin(angle), 0.0};
    // Keep the ramp and its gateways out of the bands next to the side faces
    // where horizontal-seam gateways live, when the region is wide enough.
    const double side_band = d.corridor + 2.0 * d.gateway_radius + 2.0 * d.radius + d.gap;
    auto pick = [&](double extent) {
        const double inner = reach + side_band;
        const double lo = 2.0 * inner < extent ? inner : reach;
        return lo + (extent - 2.0 * lo) * rng.uniform();
    };
    const double sx = pick(s.size.x);
    const double sy = pick(s.size.y);
    seam.position = {sx, sy, positive ? s.size.z : 0.0};
    Vec3 facing;
    if (!stable_unit(positive ? dir : Vec3{-dir.x, -dir.y, 0.0}, facing))
        throw Error(ErrorCode::constraint, "vertical seam facing has no storage-stable unit form");
    seam.facing = facing;
    seam.port = positive ? Vec3{sx - run * dir.x, sy - run * dir.y, s.size.z - delta}
                         : Vec3{sx + run * dir.x, sy + run * dir.y, delta};
    seam.slope_degrees = segment_slope_degrees(seam.port, seam.position);
    if (seam.slope_degrees > s.maximum_slope_degrees - slope_margin_degrees)
        throw Error(ErrorCode::constraint, "vertical face ramp exceeds the slope limit");
    return seam;
}

} // namespace fbs::caves::detail
