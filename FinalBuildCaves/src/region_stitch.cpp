// Cross-region stitching.
//
// Each region owns a gateway chamber with a real seam port per enabled face.
// stitch() re-derives the canonical seam from the shared settings recorded on
// both gateways, checks every shared value bit-for-bit and returns one ordinary
// tunnel from the first region's seam port, through the common face point, to
// the second region's seam port, expressed in the first region's local frame.
// The tunnel id is canonical (independent of argument order).

#include "region_internal.hpp"

namespace fbs::caves {
namespace detail {
namespace {

[[noreturn]] void reject(const std::string& message) { throw Error(ErrorCode::incompatible_boundary, message); }

struct Side {
    const Boundary* boundary = nullptr;
    const Port* port = nullptr;
    const Cavern* gateway = nullptr;
};

Side find_side(const Region& g, Face face, const char* which) {
    Side s;
    for (const Boundary& b : g.boundaries)
        if (b.face == face) s.boundary = &b;
    if (!s.boundary) reject(std::string(which) + " region has no boundary on face " + face_label(face));
    const std::string pk = uuid_key(s.boundary->port_id);
    for (const Port& p : g.ports)
        if (uuid_key(p.id) == pk) s.port = &p;
    if (!s.port) reject(std::string(which) + " region boundary references a missing port");
    const std::string ck = uuid_key(s.port->cavern_id);
    for (const Cavern& c : g.caverns)
        if (uuid_key(c.id) == ck) s.gateway = &c;
    if (!s.gateway || s.gateway->authoring != Authoring::generated || !is_generated(s.gateway->generation))
        reject(std::string(which) + " region boundary port does not belong to a generated gateway");
    if (s.port->state == PortState::sealed) reject(std::string(which) + " region seam port is sealed");
    for (const Tunnel& t : g.tunnels)
        if (uuid_key(t.start_port_id) == pk || uuid_key(t.end_port_id) == pk)
            reject(std::string(which) + " region seam port is already used inside the region");
    return s;
}

Shared side_settings(const Region& g, const Side& s, const char* which) {
    Shared shared;
    if (!shared_from_parameters(s.gateway->generation.parameters, shared))
        reject(std::string(which) + " region gateway lacks its generation settings");
    if (settings_id(shared) != g.settings_id) reject(std::string(which) + " region settings do not match its settings_id");
    if (shared.world_seed != g.world_seed || !same(shared.size, g.size))
        reject(std::string(which) + " region seed or size does not match its settings");
    if (stable_id(region_identity(shared, g.key)) != g.id) reject(std::string(which) + " region id does not match its key");
    return shared;
}

} // namespace
} // namespace detail

std::vector<Tunnel> stitch(const Region& first, const Region& second) {
    using namespace detail;
    if (first.version != schema_version || second.version != schema_version || first.algorithm != algorithm_version ||
        second.algorithm != algorithm_version)
        throw Error(ErrorCode::incompatible_version, "stitch requires regions from " + std::string(algorithm_version));
    if (first.world_seed != second.world_seed) reject("regions were generated with different world seeds");
    if (first.settings_id != second.settings_id) reject("regions were generated with different settings");
    if (!same(first.size, second.size)) reject("regions have different sizes");

    Face face = Face::negative_x;
    bool adjacent = false;
    for (int f = 0; f < 6 && !adjacent; ++f) {
        RegionKey n;
        if (neighbor_key(first.key, static_cast<Face>(f), n) && n.x == second.key.x && n.y == second.key.y &&
            n.z == second.key.z) {
            face = static_cast<Face>(f);
            adjacent = true;
        }
    }
    if (!adjacent) reject("regions are not face-adjacent");
    const Face other = opposite(face);
    const Side a = find_side(first, face, "first");
    const Side b = find_side(second, other, "second");
    const Shared shared = side_settings(first, a, "first");
    const Shared shared_b = side_settings(second, b, "second");
    if (shared_parameters(shared) != shared_parameters(shared_b)) reject("regions record different shared settings");

    const Derived derived = derive(shared);
    const Seam sa = derive_seam(shared, derived, first.key, face);
    const Seam sb = derive_seam(shared, derived, second.key, other);
    if (sa.key != sb.key || sa.width != sb.width || sa.height != sb.height) reject("canonical seams disagree");
    const Boundary& ba = *a.boundary;
    const Boundary& bb = *b.boundary;
    if (ba.key != sa.key || bb.key != sb.key) reject("boundary keys do not match the canonical seam");
    if (ba.width != sa.width || ba.height != sa.height || bb.width != sb.width || bb.height != sb.height)
        reject("boundary dimensions do not match the canonical seam");
    if (!same(ba.position, sa.position) || !same(bb.position, sb.position))
        reject("boundary positions do not match the canonical seam");
    if (!same(a.port->position, sa.port) || !same(b.port->position, sb.port) || a.port->width != sa.width ||
        a.port->height != sa.height || b.port->width != sb.width || b.port->height != sb.height)
        reject("seam ports do not match the canonical seam");

    // Translate the second region into the first region's frame: exactly one
    // axis differs by one region size.
    const int axis = face_axis(face);
    const double extent = axis == 0 ? first.size.x : axis == 1 ? first.size.y : first.size.z;
    const double shift = face_positive(face) ? extent : -extent;
    auto to_first = [&](Vec3 v) {
        (axis == 0 ? v.x : axis == 1 ? v.y : v.z) += shift;
        return v;
    };
    if (!same(to_first(bb.position), ba.position)) reject("shared seam point differs between the regions");

    const Vec3 p0 = a.port->position, p1 = ba.position, p2 = to_first(b.port->position);
    const double slope = std::max(segment_slope_degrees(p0, p1), segment_slope_degrees(p1, p2));
    if (slope > shared.maximum_slope_degrees - slope_margin_degrees)
        throw Error(ErrorCode::constraint, "seam passage would exceed the slope limit");

    const bool first_lower = face_positive(face);
    const std::string lower_port = uuid_key(first_lower ? a.port->id : b.port->id);
    const std::string upper_port = uuid_key(first_lower ? b.port->id : a.port->id);
    const RegionKey& lower = first_lower ? first.key : second.key;
    const std::string identity = "fbs-caves/1/stitch/" + uuid_key(sa.key) + "/" + lower_port + "/" + upper_port;
    Mixer128 mixer("fbs.caves.stitch.v1");
    mixer.str(identity);
    const Hash128 hash = mixer.finish();

    Tunnel t;
    t.id = stable_id(identity);
    t.start_port_id = a.port->id;
    t.end_port_id = b.port->id;
    static const char axes[] = {'X', 'Y', 'Z'};
    t.code = std::string("SEAM-") + axes[axis] + "-" + key_label(lower);
    t.name = "Seam passage " + std::string(1, axes[axis]) + " (" + key_label(lower) + ")";
    t.centerline = {{p0, sa.width, sa.height}, {p1, sa.width, sa.height}, {p2, sa.width, sa.height}};
    t.maximum_slope_degrees = shared.maximum_slope_degrees;
    t.minimum_clearance = sa.height;
    t.traversal = slope <= walkable_slope_degrees ? Traversal::walkable : Traversal::difficult;
    t.locked = false;
    t.generation.seed = static_cast<std::int32_t>(static_cast<std::uint32_t>(hash.hi >> 33));
    for (const auto& [k, v] : shared_parameters(shared)) t.generation.parameters[k] = v;
    t.generation.parameters[generator_key] = algorithm_version;
    t.generation.parameters["fbs.role"] = "seam";
    t.generation.parameters["fbs.seam_key"] = sa.key;
    t.generation.parameters["fbs.frame_region_id"] = first.id;
    t.generation.parameters["fbs.lower_region_key"] = key_label(lower);
    t.tags = {"generated", "seam"};
    return {t};
}

} // namespace fbs::caves
