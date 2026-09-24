// Input validation and Region::validate.
//
// validate() is a superset of the MOOCoW CaveDatasetValidator rules for the
// records of one region (identity, codes, references, polygons, ports,
// tunnels, transitions, connectivity) plus generator guarantees: exact
// port/profile endpoints for generated passages, floor/port elevation, slope of
// floors and centrelines, non-overlap of generated passages and chambers,
// containment in the region and exact re-derivation of every seam boundary.

#include "region_internal.hpp"
#include "utf8.hpp"

#include <algorithm>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>

namespace fbs::caves::detail {
namespace {

[[noreturn]] void fail(ErrorCode code, const std::string& message) { throw Error(code, message); }

bool blank(const std::string& s) {
    for (char c : s)
        if (!(c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\f' || c == '\v')) return false;
    return true;
}

std::string code_key(const std::string& s) {
    std::size_t b = 0, e = s.size();
    while (b < e && (s[b] == ' ' || s[b] == '\t' || s[b] == '\n' || s[b] == '\r')) ++b;
    while (e > b && (s[e - 1] == ' ' || s[e - 1] == '\t' || s[e - 1] == '\n' || s[e - 1] == '\r')) --e;
    std::string out = s.substr(b, e - b);
    for (char& c : out)
        if (c >= 'A' && c <= 'Z') c = static_cast<char>(c - 'A' + 'a');
    return out;
}

// First / last Unicode scalar of UTF-8 text (malformed bytes are returned raw).
std::uint32_t first_scalar(const std::string& s) {
    const auto* p = reinterpret_cast<const unsigned char*>(s.data());
    const std::size_t n = s.size();
    if (n == 0) return 0;
    if (p[0] < 0x80) return p[0];
    if ((p[0] & 0xe0) == 0xc0 && n >= 2) return ((p[0] & 0x1fu) << 6) | (p[1] & 0x3fu);
    if ((p[0] & 0xf0) == 0xe0 && n >= 3) return ((p[0] & 0x0fu) << 12) | ((p[1] & 0x3fu) << 6) | (p[2] & 0x3fu);
    return p[0];
}
std::uint32_t last_scalar(const std::string& s) {
    std::size_t i = s.size();
    if (i == 0) return 0;
    std::size_t start = i - 1;
    while (start > 0 && i - start < 4 && (static_cast<unsigned char>(s[start]) & 0xc0) == 0x80) --start;
    return first_scalar(s.substr(start));
}
// .NET char.IsWhiteSpace for the scalars MOOCoW's Trim() would remove.
bool dotnet_space(std::uint32_t c) {
    return (c >= 0x09 && c <= 0x0d) || c == 0x20 || c == 0x85 || c == 0xa0 || c == 0x1680 ||
           (c >= 0x2000 && c <= 0x200a) || c == 0x2028 || c == 0x2029 || c == 0x202f || c == 0x205f || c == 0x3000;
}

struct Checker {
    ErrorCode code;
    Budget* budget = nullptr;
    std::size_t vertices = 0;
    std::size_t vertex_cap = max_total_vertices;

    void work(std::uint64_t units) {
        if (budget) budget->charge(units);
    }
    // Text MOOCoW would store unchanged. Plain CSV cells (codes, names, material,
    // tags) must not contain line breaks: its CSV reader splits lines before quote
    // parsing. Generation parameters travel JSON-escaped, so only their keys are
    // held to the trimming rule. Required identifiers/names must not have
    // surrounding whitespace (MOOCoW trims with char.IsWhiteSpace).
    void text(const std::string& s, const std::string& what, bool required, bool plain = true) {
        if (s.size() > max_string_bytes) fail(code, what + " is longer than " + std::to_string(max_string_bytes) + " bytes");
        if (!valid_utf8(s)) fail(code, what + " is not valid UTF-8");
        if (plain && s.find_first_of("\r\n") != std::string::npos) fail(code, what + " must not contain line breaks");
        if (!required) return;
        if (blank(s)) fail(code, what + " must not be blank");
        if (dotnet_space(first_scalar(s)) || dotnet_space(last_scalar(s)))
            fail(code, what + " must not start or end with whitespace");
    }
    void id(const std::string& s, const std::string& what) {
        if (!is_uuid(s)) fail(code, what + " '" + s.substr(0, 64) + "' is not a non-nil 8-4-4-4-12 UUID");
    }
    void number(double v, const std::string& what) {
        if (!std::isfinite(v)) fail(code, what + " must be finite");
    }
    void positive(double v, const std::string& what) {
        number(v, what);
        if (v <= 0) fail(code, what + " must be greater than zero");
    }
    void nonnegative(double v, const std::string& what) {
        number(v, what);
        if (v < 0) fail(code, what + " must not be negative");
    }
    void vec(Vec3 v, const std::string& what) {
        if (!finite(v)) fail(code, what + " must be finite");
    }
    void tags(const std::vector<std::string>& t, const std::string& what) {
        if (t.size() > max_tags) fail(code, what + " has too many tags");
        std::set<std::string> seen;
        for (const auto& s : t) {
            text(s, what + " tag", true);
            // MOOCoW de-duplicates tags case-insensitively; a duplicate would not survive storage.
            if (!seen.insert(code_key(s)).second) fail(code, what + " has duplicate tags (case-insensitive)");
        }
        work(t.size());
    }
    void generation(const Generation& g, const std::string& what) {
        if (g.parameters.size() > max_parameters) fail(code, what + " has too many generation parameters");
        for (const auto& [k, v] : g.parameters) {
            text(k, what + " generation parameter key", true, false);
            // A root "parameters" key is read by MOOCoW as a legacy wrapper.
            if (k == "parameters") fail(code, what + " generation parameter key 'parameters' is reserved by MOOCoW");
            text(v, what + " generation parameter value", false, false);
        }
        work(g.parameters.size());
    }
    void bounds(const Bounds& b, const std::string& what) {
        vec(b.minimum, what + " minimum");
        vec(b.maximum, what + " maximum");
        if (b.minimum.x > b.maximum.x || b.minimum.y > b.maximum.y || b.minimum.z > b.maximum.z)
            fail(code, what + " minimum exceeds maximum");
        if (!std::isfinite(b.maximum.x - b.minimum.x) || !std::isfinite(b.maximum.y - b.minimum.y) ||
            !std::isfinite(b.maximum.z - b.minimum.z))
            fail(code, what + " extent overflows");
    }
    void polygon(const Polygon& p, const std::string& what) {
        if (p.size() > max_polygon_vertices) fail(code, what + " has too many vertices");
        work(1 + p.size() * p.size() / 16);
        vertices += p.size();
        if (vertices > vertex_cap) fail(code, "too many polygon/path vertices in total");
        for (const auto& v : p)
            if (!finite(v)) fail(code, what + " has a non-finite vertex");
        if (!polygon_is_simple(p)) fail(code, what + " is not a simple polygon (duplicates, crossings or degenerate)");
    }
    void path(const std::vector<Vec2>& p, const std::string& what) {
        if (p.size() > max_polygon_vertices) fail(code, what + " has too many points");
        vertices += p.size();
        if (vertices > vertex_cap) fail(code, "too many polygon/path vertices in total");
        if (p.size() < 2) fail(code, what + " needs at least two points");
        for (const auto& v : p)
            if (!finite(v)) fail(code, what + " has a non-finite point");
        if (has_consecutive_duplicates(p)) fail(code, what + " has consecutive duplicate points");
    }
};

bool inside_xy(const Bounds& b, Vec2 p) {
    return p.x >= b.minimum.x && p.x <= b.maximum.x && p.y >= b.minimum.y && p.y <= b.maximum.y;
}
bool inside(const Bounds& b, Vec3 p) {
    return p.x >= b.minimum.x && p.x <= b.maximum.x && p.y >= b.minimum.y && p.y <= b.maximum.y &&
           p.z >= b.minimum.z && p.z <= b.maximum.z;
}
bool has_tag(const std::vector<std::string>& tags, const char* tag) {
    for (const auto& t : tags)
        if (code_key(t) == tag) return true;
    return false;
}

// Record-level rules shared by requests (authored records) and regions.
struct Records {
    const std::vector<Cavern>& caverns;
    const std::vector<Floor>& floors;
    const std::vector<Cliff>& cliffs;
    const std::vector<Port>& ports;
    const std::vector<Tunnel>& tunnels;
};

struct Index {
    std::unordered_map<std::string, std::size_t> cavern, floor, port;
};

Index check_records(const Records& r, Checker& c, std::unordered_set<std::string>& ids) {
    Index ix;
    auto unique_id = [&](const std::string& id, const std::string& what) {
        c.id(id, what + " id");
        if (!ids.insert(uuid_key(id)).second) fail(c.code, what + " id " + id + " is duplicated");
    };
    // MOOCoW code scopes: caverns and tunnels share the map scope; floors, cliffs
    // and ports share their cavern's scope (case-insensitive, trimmed).
    std::set<std::string> map_codes;
    std::set<std::pair<std::string, std::string>> cavern_codes;

    for (std::size_t i = 0; i < r.caverns.size(); ++i) {
        const Cavern& v = r.caverns[i];
        const std::string what = "cavern " + v.id;
        unique_id(v.id, "cavern");
        c.text(v.code, what + " code", true);
        c.text(v.name, what + " name", true);
        if (!map_codes.insert(code_key(v.code)).second) fail(c.code, "cavern code '" + v.code + "' is duplicated in the map scope");
        c.work(4);
        c.vec(v.center, what + " center");
        c.bounds(v.bounds, what + " bounds");
        if (!inside(v.bounds, v.center)) fail(c.code, what + " center is outside its bounds");
        c.nonnegative(v.reserved_clearance, what + " reserved clearance");
        if (!v.boundary.empty()) {
            c.polygon(v.boundary, what + " boundary");
            for (const auto& p : v.boundary)
                if (!inside_xy(v.bounds, p)) fail(c.code, what + " boundary is outside its bounds");
        } else if (v.authoring == Authoring::exact) {
            fail(c.code, what + " is exact but has no boundary");
        }
        c.generation(v.generation, what);
        c.tags(v.tags, what);
        ix.cavern.emplace(uuid_key(v.id), i);
    }
    for (std::size_t i = 0; i < r.floors.size(); ++i) {
        const Floor& f = r.floors[i];
        const std::string what = "floor " + f.id;
        unique_id(f.id, "floor");
        c.work(4);
        c.text(f.code, what + " code", true);
        c.text(f.name, what + " name", true);
        c.text(f.material, what + " material", true);
        const auto it = ix.cavern.find(uuid_key(f.cavern_id));
        if (it == ix.cavern.end()) fail(c.code, what + " references a missing cavern");
        if (!cavern_codes.insert({uuid_key(f.cavern_id), code_key(f.code)}).second)
            fail(c.code, what + " code is duplicated in its cavern");
        c.polygon(f.boundary, what + " boundary");
        for (const auto& p : f.boundary)
            if (!inside_xy(r.caverns[it->second].bounds, p)) fail(c.code, what + " is outside its cavern bounds");
        c.number(f.base_z, what + " base elevation");
        c.number(f.slope_x, what + " slope x");
        c.number(f.slope_y, what + " slope y");
        c.nonnegative(f.variation_amplitude, what + " variation amplitude");
        c.tags(f.tags, what);
        ix.floor.emplace(uuid_key(f.id), i);
    }
    for (const Cliff& k : r.cliffs) {
        const std::string what = "cliff " + k.id;
        unique_id(k.id, "cliff");
        c.work(4 + k.boundary.size());
        c.text(k.code, what + " code", true);
        if (!ix.cavern.count(uuid_key(k.cavern_id))) fail(c.code, what + " references a missing cavern");
        if (!cavern_codes.insert({uuid_key(k.cavern_id), code_key(k.code)}).second)
            fail(c.code, what + " code is duplicated in its cavern");
        const auto from = ix.floor.find(uuid_key(k.from_floor_id));
        if (from == ix.floor.end()) fail(c.code, what + " references a missing from-floor");
        if (uuid_key(r.floors[from->second].cavern_id) != uuid_key(k.cavern_id))
            fail(c.code, what + " from-floor belongs to another cavern");
        if (!k.to_floor_id.empty()) {
            const auto to = ix.floor.find(uuid_key(k.to_floor_id));
            if (to == ix.floor.end()) fail(c.code, what + " references a missing to-floor");
            if (uuid_key(r.floors[to->second].cavern_id) != uuid_key(k.cavern_id))
                fail(c.code, what + " to-floor belongs to another cavern");
            if (to->second == from->second) fail(c.code, what + " connects a floor to itself");
        }
        c.path(k.boundary, what + " boundary");
        c.nonnegative(k.height, what + " height");
        if (k.transition == Transition::one_way_drop && k.to_floor_id.empty())
            fail(c.code, what + " one-way drop needs a destination floor");
        const bool needs_height = k.transition == Transition::slope || k.transition == Transition::ramp ||
                                  k.transition == Transition::stairs || k.transition == Transition::climbable_cliff ||
                                  k.transition == Transition::impassable_cliff || k.transition == Transition::one_way_drop;
        if (needs_height && k.height <= geometry_tolerance) fail(c.code, what + " transition needs a positive height");
    }
    for (std::size_t i = 0; i < r.ports.size(); ++i) {
        const Port& p = r.ports[i];
        const std::string what = "port " + p.id;
        unique_id(p.id, "port");
        c.work(4);
        c.text(p.code, what + " code", true);
        c.text(p.name, what + " name", true);
        const auto cav = ix.cavern.find(uuid_key(p.cavern_id));
        if (cav == ix.cavern.end()) fail(c.code, what + " references a missing cavern");
        if (!cavern_codes.insert({uuid_key(p.cavern_id), code_key(p.code)}).second)
            fail(c.code, what + " code is duplicated in its cavern");
        c.vec(p.position, what + " position");
        c.vec(p.facing, what + " facing");
        if (length(p.facing) <= 1e-12) fail(c.code, what + " facing must be non-zero");
        c.positive(p.width, what + " width");
        c.positive(p.height, what + " height");
        if (!inside(r.caverns[cav->second].bounds, p.position)) fail(c.code, what + " is outside its cavern bounds");
        if (!p.floor_id.empty()) {
            const auto fl = ix.floor.find(uuid_key(p.floor_id));
            if (fl == ix.floor.end()) fail(c.code, what + " references a missing floor");
            const Floor& f = r.floors[fl->second];
            if (uuid_key(f.cavern_id) != uuid_key(p.cavern_id)) fail(c.code, what + " floor belongs to another cavern");
            if (!polygon_contains(f.boundary, {p.position.x, p.position.y}))
                fail(c.code, what + " is outside its floor region");
        }
        ix.port.emplace(uuid_key(p.id), i);
    }
    for (const Tunnel& t : r.tunnels) {
        const std::string what = "tunnel " + t.id;
        unique_id(t.id, "tunnel");
        c.text(t.code, what + " code", true);
        c.text(t.name, what + " name", true);
        if (!map_codes.insert(code_key(t.code)).second) fail(c.code, "tunnel code '" + t.code + "' is duplicated in the map scope");
        c.work(4 + t.centerline.size());
        const auto a = ix.port.find(uuid_key(t.start_port_id));
        const auto b = ix.port.find(uuid_key(t.end_port_id));
        if (a == ix.port.end() || b == ix.port.end()) fail(c.code, what + " references a missing port");
        if (a->second == b->second) fail(c.code, what + " must connect two different ports");
        if (t.centerline.size() < 2) fail(c.code, what + " needs at least two profiles");
        if (t.centerline.size() > max_tunnel_profiles) fail(c.code, what + " has too many profiles");
        c.vertices += t.centerline.size();
        if (c.vertices > c.vertex_cap) fail(c.code, "too many polygon/path vertices in total");
        c.number(t.maximum_slope_degrees, what + " maximum slope");
        if (t.maximum_slope_degrees < 0 || t.maximum_slope_degrees >= 90)
            fail(c.code, what + " maximum slope must be in [0, 90) degrees");
        c.positive(t.minimum_clearance, what + " minimum clearance");
        for (const Profile& pr : t.centerline) {
            c.vec(pr.center, what + " profile");
            c.positive(pr.width, what + " profile width");
            c.positive(pr.height, what + " profile height");
            if (pr.height < t.minimum_clearance) fail(c.code, what + " profile height is below its minimum clearance");
        }
        for (std::size_t i = 1; i < t.centerline.size(); ++i)
            if (segment_slope_degrees(t.centerline[i - 1].center, t.centerline[i].center) >
                t.maximum_slope_degrees + geometry_tolerance)
                fail(c.code, what + " exceeds its maximum slope");
        const Port& sp = r.ports[a->second];
        const Port& ep = r.ports[b->second];
        if (distance(t.centerline.front().center, sp.position) > std::max(sp.width, sp.height))
            fail(c.code, what + " misses its start port");
        if (distance(t.centerline.back().center, ep.position) > std::max(ep.width, ep.height))
            fail(c.code, what + " misses its end port");
        if (sp.state == PortState::sealed || ep.state == PortState::sealed) fail(c.code, what + " uses a sealed port");
        c.generation(t.generation, what);
        c.tags(t.tags, what);
    }
    return ix;
}

// MOOCoW floor connectivity: floors of one cavern must be connected through
// non-blocking transitions unless the isolated component is tagged "isolated".
void check_floor_connectivity(const Records& r, const Index& ix, ErrorCode code) {
    std::unordered_map<std::string, std::vector<std::size_t>> by_cavern;
    for (std::size_t i = 0; i < r.floors.size(); ++i) by_cavern[uuid_key(r.floors[i].cavern_id)].push_back(i);
    std::unordered_map<std::size_t, std::vector<std::size_t>> links;
    for (const Cliff& k : r.cliffs) {
        if (k.to_floor_id.empty() || k.transition == Transition::unconnected ||
            k.transition == Transition::impassable_cliff)
            continue;
        const std::size_t a = ix.floor.at(uuid_key(k.from_floor_id)), b = ix.floor.at(uuid_key(k.to_floor_id));
        links[a].push_back(b);
        links[b].push_back(a);
    }
    for (auto& [cavern, list] : by_cavern) {
        if (list.size() < 2) continue;
        std::set<std::size_t> remaining(list.begin(), list.end());
        while (!remaining.empty()) {
            std::vector<std::size_t> component, stack{*remaining.begin()};
            remaining.erase(remaining.begin());
            while (!stack.empty()) {
                const std::size_t u = stack.back();
                stack.pop_back();
                component.push_back(u);
                for (std::size_t v : links[u])
                    if (remaining.erase(v)) stack.push_back(v);
            }
            if (component.size() >= list.size()) continue;
            for (std::size_t f : component)
                if (!has_tag(r.floors[f].tags, "isolated"))
                    fail(code, "floor " + r.floors[f].id + " is in an isolated component without an 'isolated' tag");
        }
    }
}

double tunnel_radius(const Tunnel& t) {
    double w = 0, h = 0;
    for (const auto& p : t.centerline) {
        w = std::max(w, p.width);
        h = std::max(h, p.height);
    }
    return section_radius(w, h);
}

Bounds tunnel_box(const Tunnel& t, double r) {
    Bounds b{t.centerline.front().center, t.centerline.front().center};
    for (const auto& p : t.centerline) {
        b.minimum = {std::min(b.minimum.x, p.center.x), std::min(b.minimum.y, p.center.y), std::min(b.minimum.z, p.center.z)};
        b.maximum = {std::max(b.maximum.x, p.center.x), std::max(b.maximum.y, p.center.y), std::max(b.maximum.z, p.center.z)};
    }
    return inflate(b, r);
}

bool segment_near_box(Vec3 a, Vec3 b, const Bounds& box, double threshold) {
    Bounds s{{std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z)},
             {std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z)}};
    if (!boxes_overlap(inflate(s, threshold), box)) return false;
    return segment_box_distance(a, b, box) < threshold;
}

} // namespace

void validate_request(const Request& q, const Options& o, Budget* budget) {
    Checker c{ErrorCode::invalid_input, budget};
    if (q.version != schema_version)
        fail(ErrorCode::incompatible_version, "request schema version " + std::to_string(q.version) + " is not supported");
    if (!finite(q.size)) fail(c.code, "region size must be finite");
    for (double v : {q.size.x, q.size.y, q.size.z})
        if (v < min_region_extent || v > max_region_extent)
            fail(c.code, "region size must be between " + format_double(min_region_extent) + " and " +
                             format_double(max_region_extent) + " m per axis");
    if (q.face_mask > 63) fail(c.code, "face_mask has bits outside the six region faces");
    if (q.chamber_count > max_chambers) fail(c.code, "chamber_count exceeds " + std::to_string(max_chambers));
    if (q.loop_count > max_loops) fail(c.code, "loop_count exceeds " + std::to_string(max_loops));
    c.positive(q.chamber_radius_min, "chamber_radius_min");
    c.positive(q.chamber_radius_max, "chamber_radius_max");
    if (q.chamber_radius_min > q.chamber_radius_max) fail(c.code, "chamber_radius_min exceeds chamber_radius_max");
    c.positive(q.chamber_height, "chamber_height");
    c.positive(q.tunnel_width, "tunnel_width");
    c.positive(q.tunnel_height, "tunnel_height");
    for (double v : {q.chamber_radius_max, q.chamber_height, q.tunnel_width, q.tunnel_height})
        if (v > max_feature_size) fail(c.code, "feature dimensions must not exceed " + format_double(max_feature_size) + " m");
    c.number(q.maximum_slope_degrees, "maximum_slope_degrees");
    if (q.maximum_slope_degrees <= 0 || q.maximum_slope_degrees >= 90)
        fail(c.code, "maximum_slope_degrees must be in (0, 90)");
    c.positive(q.sampling_radius, "sampling_radius");
    if (q.sampling_radius > 1) fail(c.code, "sampling_radius must not exceed 1");
    if (q.neighbors < 3 || q.neighbors > max_neighbors)
        fail(c.code, "neighbors must be between 3 and " + std::to_string(max_neighbors));
    if (q.reserved.size() > max_reserved_boxes) fail(c.code, "too many reserved boxes");
    for (const Bounds& b : q.reserved) c.bounds(b, "reserved box");
    for (std::size_t n : {q.authored_caverns.size(), q.authored_floors.size(), q.authored_cliffs.size(),
                          q.authored_ports.size(), q.authored_tunnels.size()})
        if (n > max_authored_records) fail(c.code, "too many authored records of one kind");
    if (o.work_limit == 0) fail(c.code, "options.work_limit must be positive");
    if (o.max_points < 16 || o.max_edges < 16) fail(c.code, "options.max_points and max_edges must be at least 16");

    std::unordered_set<std::string> ids;
    const Records r{q.authored_caverns, q.authored_floors, q.authored_cliffs, q.authored_ports, q.authored_tunnels};
    const Index ix = check_records(r, c, ids);
    check_floor_connectivity(r, ix, c.code);

    // Authored geometry must lie inside the region (local frame); cross-region
    // content is expressed through seams, never by records spilling over a face.
    const Bounds region{{0, 0, 0}, q.size};
    for (const Cavern& v : q.authored_caverns)
        if (!inside(region, v.bounds.minimum) || !inside(region, v.bounds.maximum))
            fail(c.code, "authored cavern " + v.id + " bounds are outside the region");
    for (const Tunnel& t : q.authored_tunnels)
        for (const Profile& p : t.centerline)
            if (!inside(region, p.center)) fail(c.code, "authored tunnel " + t.id + " leaves the region");
}

void check_region(const Region& g, Budget* budget) {
    const ErrorCode code = ErrorCode::invalid_input;
    if (g.version != schema_version) fail(ErrorCode::incompatible_version, "region schema version is not supported");
    if (g.algorithm != algorithm_version) fail(ErrorCode::incompatible_version, "region algorithm is not " + std::string(algorithm_version));
    if (!is_uuid(g.id) || !is_uuid(g.settings_id)) fail(code, "region id and settings_id must be UUIDs");
    if (!finite(g.size) || !finite(g.origin)) fail(code, "region size and origin must be finite");
    for (double v : {g.size.x, g.size.y, g.size.z})
        if (v < min_region_extent || v > max_region_extent) fail(code, "region size is out of range");
    for (std::size_t n : {g.caverns.size(), g.floors.size(), g.cliffs.size(), g.ports.size(), g.tunnels.size()})
        if (n > max_region_records) fail(code, "region has too many records");
    if (g.boundaries.size() > 6) fail(code, "region has more than six boundaries");

    Checker c{code, budget};
    c.vertex_cap = 4 * max_total_vertices;
    std::unordered_set<std::string> ids{uuid_key(g.id)};
    const Records r{g.caverns, g.floors, g.cliffs, g.ports, g.tunnels};
    const Index ix = check_records(r, c, ids);
    check_floor_connectivity(r, ix, code);

    // Shared settings recorded on generated caverns must reproduce settings_id and id.
    std::vector<bool> generated_cavern(g.caverns.size(), false);
    bool have_shared = false;
    Shared shared;
    for (std::size_t i = 0; i < g.caverns.size(); ++i) {
        const Cavern& v = g.caverns[i];
        if (v.authoring != Authoring::generated || !is_generated(v.generation)) continue;
        generated_cavern[i] = true;
        Shared s;
        if (!shared_from_parameters(v.generation.parameters, s)) fail(code, "generated cavern " + v.id + " lacks shared settings");
        if (settings_id(s) != g.settings_id) fail(code, "generated cavern " + v.id + " settings do not match settings_id");
        if (!have_shared) {
            shared = s;
            have_shared = true;
        }
    }
    if (have_shared) {
        if (shared.world_seed != g.world_seed) fail(code, "region world_seed does not match its settings");
        if (!same(shared.size, g.size)) fail(code, "region size does not match its settings");
        if (stable_id(region_identity(shared, g.key)) != g.id) fail(code, "region id does not match its key and settings");
        const double ox = static_cast<double>(g.key.x) * g.size.x, oy = static_cast<double>(g.key.y) * g.size.y,
                     oz = static_cast<double>(g.key.z) * g.size.z;
        if (!(g.origin.x == ox && g.origin.y == oy && g.origin.z == oz)) fail(code, "region origin does not match key * size");
    }

    std::vector<bool> generated_tunnel(g.tunnels.size());
    for (std::size_t i = 0; i < g.tunnels.size(); ++i) generated_tunnel[i] = is_generated(g.tunnels[i].generation);
    auto port_cavern = [&](const std::string& port_id) {
        return ix.cavern.at(uuid_key(g.ports[ix.port.at(uuid_key(port_id))].cavern_id));
    };

    // Port usage (required / sealed), generated endpoint exactness and floor slope.
    std::unordered_map<std::string, int> used;
    for (std::size_t i = 0; i < g.tunnels.size(); ++i) {
        const Tunnel& t = g.tunnels[i];
        ++used[uuid_key(t.start_port_id)];
        ++used[uuid_key(t.end_port_id)];
        if (!generated_tunnel[i]) continue;
        const Port& sp = g.ports[ix.port.at(uuid_key(t.start_port_id))];
        const Port& ep = g.ports[ix.port.at(uuid_key(t.end_port_id))];
        const Profile& f = t.centerline.front();
        const Profile& l = t.centerline.back();
        if (!same(f.center, sp.position) || f.width != sp.width || f.height != sp.height)
            fail(code, "generated tunnel " + t.id + " does not start exactly at its port");
        if (!same(l.center, ep.position) || l.width != ep.width || l.height != ep.height)
            fail(code, "generated tunnel " + t.id + " does not end exactly at its port");
        for (std::size_t k = 1; k < t.centerline.size(); ++k) {
            const Vec3 a = t.centerline[k - 1].center, b = t.centerline[k].center;
            const Vec3 fa{a.x, a.y, a.z - 0.5 * t.centerline[k - 1].height};
            const Vec3 fb{b.x, b.y, b.z - 0.5 * t.centerline[k].height};
            if (segment_slope_degrees(fa, fb) > t.maximum_slope_degrees + geometry_tolerance)
                fail(code, "generated tunnel " + t.id + " floor exceeds its maximum slope");
            if (!inside({{0, 0, 0}, g.size}, b) || !inside({{0, 0, 0}, g.size}, a))
                fail(code, "generated tunnel " + t.id + " leaves the region");
        }
    }
    for (const Port& p : g.ports) {
        const int n = used.count(uuid_key(p.id)) ? used[uuid_key(p.id)] : 0;
        if (p.state == PortState::required && n == 0) fail(code, "required port " + p.id + " is not connected");
        if (p.state == PortState::sealed && n != 0) fail(code, "sealed port " + p.id + " is used");
        const std::size_t cv = ix.cavern.at(uuid_key(p.cavern_id));
        if (generated_cavern[cv]) {
            if (n > 1) fail(code, "generated port " + p.id + " is used by more than one tunnel");
            if (!is_stable_unit(p.facing)) fail(code, "generated port " + p.id + " facing is not a storage-stable unit vector");
            if (p.floor_id.empty()) fail(code, "generated port " + p.id + " has no floor");
            const Floor& f = g.floors[ix.floor.at(uuid_key(p.floor_id))];
            if (f.slope_x != 0 || f.slope_y != 0 || f.variation_amplitude != 0 ||
                std::fabs(p.position.z - 0.5 * p.height - f.base_z) > 1e-9 * (1.0 + std::fabs(f.base_z)))
                fail(code, "generated port " + p.id + " bottom is not on its floor");
            if (p.position.z + 0.5 * p.height > g.caverns[cv].bounds.maximum.z + 1e-9)
                fail(code, "generated port " + p.id + " is taller than its cavern");
        }
    }

    // Connectivity by explicit identities: cavern graph through tunnel ports.
    std::vector<std::size_t> parent(g.caverns.size());
    for (std::size_t i = 0; i < parent.size(); ++i) parent[i] = i;
    auto find = [&](std::size_t x) {
        while (parent[x] != x) x = parent[x] = parent[parent[x]];
        return x;
    };
    std::vector<bool> linked(g.caverns.size(), false);
    for (const Tunnel& t : g.tunnels) {
        const std::size_t a = port_cavern(t.start_port_id), b = port_cavern(t.end_port_id);
        linked[a] = linked[b] = true;
        parent[find(a)] = find(b);
    }
    std::size_t main = g.caverns.size();
    for (std::size_t i = 0; i < g.caverns.size() && main == g.caverns.size(); ++i)
        if (generated_cavern[i]) main = find(i);
    for (std::size_t i = 0; i < g.caverns.size(); ++i) {
        if (has_tag(g.caverns[i].tags, "isolated")) continue;
        if (!linked[i]) fail(code, "cavern " + g.caverns[i].id + " has no connected ports and is not tagged isolated");
        if (main != g.caverns.size() && find(i) != main)
            fail(code, "cavern " + g.caverns[i].id + " is not connected to the regional network");
    }

    // Broad phase: sweep and prune along x over cavern and tunnel boxes, so the
    // overlap checks stay near-linear for well-formed regions.
    const Bounds region_box{{0, 0, 0}, g.size};
    std::vector<double> radius(g.tunnels.size());
    std::vector<Bounds> box(g.tunnels.size());
    for (std::size_t i = 0; i < g.tunnels.size(); ++i) {
        radius[i] = tunnel_radius(g.tunnels[i]);
        box[i] = tunnel_box(g.tunnels[i], radius[i]);
    }
    struct Item { double lo, hi; std::size_t index; bool tunnel; };
    std::vector<Item> items;
    items.reserve(g.caverns.size() + g.tunnels.size());
    for (std::size_t i = 0; i < g.caverns.size(); ++i)
        items.push_back({g.caverns[i].bounds.minimum.x, g.caverns[i].bounds.maximum.x, i, false});
    for (std::size_t i = 0; i < g.tunnels.size(); ++i) items.push_back({box[i].minimum.x, box[i].maximum.x, i, true});
    std::sort(items.begin(), items.end(), [](const Item& a, const Item& b) {
        return a.lo != b.lo ? a.lo < b.lo : (a.tunnel != b.tunnel ? !a.tunnel : a.index < b.index);
    });

    // Generated chambers: inside the region, with an outline, not overlapping any other cavern.
    for (std::size_t i = 0; i < g.caverns.size(); ++i) {
        if (!generated_cavern[i]) continue;
        const Cavern& v = g.caverns[i];
        if (!inside(region_box, v.bounds.minimum) || !inside(region_box, v.bounds.maximum))
            fail(code, "generated cavern " + v.id + " is outside the region");
        if (v.boundary.empty()) fail(code, "generated cavern " + v.id + " has no boundary");
    }

    // Passages: no overlap with caverns except at their own endpoint access
    // segments, and no overlap between different passages (crossings must pass
    // at distinct elevations; nothing is implicitly joined).
    auto tunnel_cavern = [&](std::size_t ti, std::size_t cv) {
        c.work(1);
        if (!generated_tunnel[ti] && !generated_cavern[cv]) return;
        const Tunnel& t = g.tunnels[ti];
        const Bounds& cb = g.caverns[cv].bounds;
        if (!boxes_overlap(box[ti], cb)) return;
        const std::size_t sc = port_cavern(t.start_port_id), ec = port_cavern(t.end_port_id);
        const std::size_t n = t.centerline.size();
        c.work(n);
        for (std::size_t k = 1; k < n; ++k) {
            if (k == 1 && cv == sc) continue;
            if (k == n - 1 && cv == ec) continue;
            if (segment_near_box(t.centerline[k - 1].center, t.centerline[k].center, cb, radius[ti] - 1e-6))
                fail(code, "tunnel " + t.id + " passes through cavern " + g.caverns[cv].id);
        }
    };
    auto tunnel_tunnel = [&](std::size_t i, std::size_t j) {
        c.work(1);
        if (!generated_tunnel[i] && !generated_tunnel[j]) return;
        if (!boxes_overlap(box[i], box[j])) return;
        c.work(g.tunnels[i].centerline.size() * g.tunnels[j].centerline.size());
        const Tunnel& t = g.tunnels[i];
        const Tunnel& u = g.tunnels[j];
        const double limit = radius[i] + radius[j] - 1e-6;
        for (std::size_t a = 1; a < t.centerline.size(); ++a) {
            const Vec3 p = t.centerline[a - 1].center, q = t.centerline[a].center;
            for (std::size_t b = 1; b < u.centerline.size(); ++b)
                if (segment_segment_distance(p, q, u.centerline[b - 1].center, u.centerline[b].center) < limit)
                    fail(code, "tunnels " + t.id + " and " + u.id + " overlap without an explicit junction");
        }
    };
    // Only pairs involving a generated item are checked (authored/authored pairs
    // are the author's responsibility), so probe from generated items only:
    // candidates are the x-sorted items whose interval can reach the probe.
    // Work is proportional to generated items x overlapping neighbours, never
    // to authored x authored pairs, and every candidate visit is charged.
    double max_extent = 0;
    for (const Item& it : items) max_extent = std::max(max_extent, it.hi - it.lo);
    auto is_generated_item = [&](const Item& it) {
        return it.tunnel ? static_cast<bool>(generated_tunnel[it.index]) : static_cast<bool>(generated_cavern[it.index]);
    };
    for (std::size_t a = 0; a < items.size(); ++a) {
        const Item& x = items[a];
        if (!is_generated_item(x)) continue;
        const double from = x.lo - max_extent;
        auto first = std::lower_bound(items.begin(), items.end(), from,
                                      [](const Item& it, double v) { return it.lo < v; });
        for (std::size_t b = static_cast<std::size_t>(first - items.begin()); b < items.size() && items[b].lo <= x.hi; ++b) {
            c.work(1);
            if (b == a) continue;
            const Item& y = items[b];
            if (y.hi < x.lo) continue;
            if (is_generated_item(y) && b < a) continue; // generated pair: handled once, from the earlier item
            if (!x.tunnel && !y.tunnel) {
                if (boxes_overlap(g.caverns[x.index].bounds, g.caverns[y.index].bounds))
                    fail(code, "generated cavern " + g.caverns[x.index].id + " overlaps cavern " + g.caverns[y.index].id);
            } else if (x.tunnel && y.tunnel) {
                tunnel_tunnel(x.index, y.index);
            } else {
                tunnel_cavern(x.tunnel ? x.index : y.index, x.tunnel ? y.index : x.index);
            }
        }
    }

    // Boundaries: one per face, own gateway port, exact canonical seam.
    std::set<int> faces;
    std::set<std::string> keys;
    for (const Boundary& b : g.boundaries) {
        if (static_cast<std::uint32_t>(b.face) > 5) fail(code, "boundary face is invalid");
        if (!faces.insert(static_cast<int>(b.face)).second) fail(code, "boundary face is duplicated");
        if (!is_uuid(b.key) || !keys.insert(uuid_key(b.key)).second) fail(code, "boundary key is invalid or duplicated");
        const auto pi_ = ix.port.find(uuid_key(b.port_id));
        if (pi_ == ix.port.end()) fail(code, "boundary references a missing port");
        const Port& p = g.ports[pi_->second];
        if (!generated_cavern[ix.cavern.at(uuid_key(p.cavern_id))]) fail(code, "boundary port must belong to a generated gateway");
        if (!have_shared) fail(code, "boundary present without generated settings");
        const Seam s = derive_seam(shared, derive(shared), g.key, b.face);
        if (s.key != b.key || !same(s.position, b.position) || s.width != b.width || s.height != b.height)
            fail(code, "boundary on face " + std::string(face_label(b.face)) + " does not match its canonical seam");
        if (!same(p.position, s.port) || p.width != s.width || p.height != s.height)
            fail(code, "boundary port on face " + std::string(face_label(b.face)) + " does not match its canonical seam");
    }
}

} // namespace fbs::caves::detail

namespace fbs::caves {
// Public validation runs under the default Options work ceiling so that
// oversized or adversarial regions end with Error(budget) instead of hanging.
void validate(const Region& region) {
    const Options options;
    detail::Budget budget(options);
    detail::check_region(region, &budget);
}
} // namespace fbs::caves
