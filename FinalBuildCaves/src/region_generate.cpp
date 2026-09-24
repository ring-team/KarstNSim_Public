// Regional cave generation.
//
// Pipeline (all state is local to one generate() call; no globals, no global RNG):
//  1. Validate the request; copy authored records verbatim.
//  2. Static exclusions: authored cavern bounds (+ their reserved clearance),
//     reserved boxes, authored tunnel capsules.
//  3. Gateways: one small chamber per enabled face, holding the canonical seam
//     port; the straight port->seam corridor is reserved.
//  4. Chambers: seeded, non-overlapping irregular star polygons with flat floors,
//     optional raised shelves joined by explicit cliff transitions.
//  5. Support lattice: a jittered lattice with slope-limited candidate edges.
//  6. Connections: spanning network (Kruskal over chamber pairs), deliberate
//     loops, required authored ports, isolated authored components. Each
//     connection gets two real ports, horizontal access stubs and a KarstNSim
//     route over the lattice subgraph whose every edge is clear of all
//     exclusions (including previously generated passages).
//  7. Region::validate plus request-level checks (reserved boxes, clearances).
// Any failure throws; nothing partial is ever returned.

#include "karst_adapter.hpp"
#include "region_internal.hpp"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <limits>
#include <tuple>
#include <map>
#include <optional>
#include <unordered_map>

namespace fbs::caves {
namespace detail {

namespace {

constexpr double port_inset = 0.25;       // port centre inside the floor polygon edge
constexpr double keypoint_spacing = 1.0;  // support nodes kept this far from access points

std::int32_t seed32(const Hash128& h) { return static_cast<std::int32_t>(static_cast<std::uint32_t>(h.hi >> 33)); }

double norm_angle(double a) {
    a = std::fmod(a, 2.0 * pi);
    return a < 0 ? a + 2.0 * pi : a;
}

// Star-shaped polygon: strictly increasing angles around a centre.
struct Star {
    Vec2 center;
    std::vector<double> angle, radius;

    Vec2 vertex(std::size_t i) const {
        return {center.x + radius[i] * std::cos(angle[i]), center.y + radius[i] * std::sin(angle[i])};
    }
    Polygon polygon() const {
        Polygon p;
        for (std::size_t i = 0; i < angle.size(); ++i) p.push_back(vertex(i));
        return p;
    }
    // Distance from the centre to the boundary along direction psi.
    double boundary(double psi) const {
        const std::size_t n = angle.size();
        const double t = norm_angle(psi - angle[0]);
        std::size_t i = 0;
        while (i + 1 < n && angle[i + 1] - angle[0] <= t) ++i;
        const Vec2 a = vertex(i), b = vertex((i + 1) % n);
        const Vec2 d{std::cos(psi), std::sin(psi)};
        const Vec2 e{b.x - a.x, b.y - a.y};
        const Vec2 w{a.x - center.x, a.y - center.y};
        const double den = d.x * e.y - d.y * e.x;
        if (std::fabs(den) < 1e-12) return std::min(radius[i], radius[(i + 1) % n]);
        return (w.x * e.y - w.y * e.x) / den;
    }
};

Star make_star(Rng& rng, Vec2 center, double r, std::size_t n) {
    Star s;
    s.center = center;
    const double theta0 = rng.uniform() * 2.0 * pi / static_cast<double>(n);
    for (std::size_t i = 0; i < n; ++i) {
        s.angle.push_back(theta0 + 2.0 * pi * static_cast<double>(i) / static_cast<double>(n));
        s.radius.push_back(r * (0.72 + 0.28 * rng.uniform()));
    }
    return s;
}

Bounds polygon_bounds(const Polygon& p, double z0, double z1) {
    Bounds b{{p[0].x, p[0].y, z0}, {p[0].x, p[0].y, z1}};
    for (const auto& v : p) {
        b.minimum.x = std::min(b.minimum.x, v.x);
        b.minimum.y = std::min(b.minimum.y, v.y);
        b.maximum.x = std::max(b.maximum.x, v.x);
        b.maximum.y = std::max(b.maximum.y, v.y);
    }
    return b;
}

Bounds segment_bounds(Vec3 a, Vec3 b) {
    return {{std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z)},
            {std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z)}};
}

double box_box_distance(const Bounds& a, const Bounds& b) {
    const double dx = std::max({0.0, a.minimum.x - b.maximum.x, b.minimum.x - a.maximum.x});
    const double dy = std::max({0.0, a.minimum.y - b.maximum.y, b.minimum.y - a.maximum.y});
    const double dz = std::max({0.0, a.minimum.z - b.maximum.z, b.minimum.z - a.maximum.z});
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// XY bucket grid over exclusion bounds. Very large exclusions go to a list that
// every query scans, so registration never exceeds a fixed number of cells.
class Grid {
public:
    void init(const Bounds& area, double cell) {
        x0_ = area.minimum.x;
        y0_ = area.minimum.y;
        cell_ = cell;
        nx_ = std::max(1, static_cast<int>(std::ceil((area.maximum.x - area.minimum.x) / cell)));
        ny_ = std::max(1, static_cast<int>(std::ceil((area.maximum.y - area.minimum.y) / cell)));
        cells_.assign(static_cast<std::size_t>(nx_) * static_cast<std::size_t>(ny_), {});
    }
    void add(std::uint32_t index, const Bounds& b) {
        int i0, i1, j0, j1;
        range(b, i0, i1, j0, j1);
        if (static_cast<long long>(i1 - i0 + 1) * (j1 - j0 + 1) > 1024) {
            large_.push_back(index);
        } else {
            for (int j = j0; j <= j1; ++j)
                for (int i = i0; i <= i1; ++i) cells_[static_cast<std::size_t>(j) * nx_ + i].push_back(index);
        }
        if (stamp_.size() <= index) stamp_.resize(index + 1, 0);
    }
    template <class F> bool all(const Bounds& q, F f) {
        ++epoch_;
        for (std::uint32_t e : large_)
            if (!f(e)) return false;
        int i0, i1, j0, j1;
        range(q, i0, i1, j0, j1);
        for (int j = j0; j <= j1; ++j)
            for (int i = i0; i <= i1; ++i)
                for (std::uint32_t e : cells_[static_cast<std::size_t>(j) * nx_ + i]) {
                    if (stamp_[e] == epoch_) continue;
                    stamp_[e] = epoch_;
                    if (!f(e)) return false;
                }
        return true;
    }

private:
    void range(const Bounds& b, int& i0, int& i1, int& j0, int& j1) const {
        auto clampi = [](double v, int n) {
            if (!(v > 0)) return 0;
            if (v >= n - 1) return n - 1;
            return static_cast<int>(v);
        };
        i0 = clampi((b.minimum.x - x0_) / cell_, nx_);
        i1 = clampi((b.maximum.x - x0_) / cell_, nx_);
        j0 = clampi((b.minimum.y - y0_) / cell_, ny_);
        j1 = clampi((b.maximum.y - y0_) / cell_, ny_);
    }
    double x0_ = 0, y0_ = 0, cell_ = 1;
    int nx_ = 1, ny_ = 1;
    std::vector<std::vector<std::uint32_t>> cells_;
    std::vector<std::uint32_t> large_;
    std::vector<std::uint32_t> stamp_;
    std::uint32_t epoch_ = 0;
};

struct Host {
    std::size_t cavern = 0;
    std::string identity, key, label, name;
    Star star;
    double floor_z = 0, ceiling_z = 0;
    bool gateway = false, shelf = false;
    double shelf_height = 0, cut_angle = 0;
    std::size_t main_floor = 0, shelf_floor = 0;  // indices into region floors
    std::vector<Vec3> ports;
    std::vector<double> port_radius;
    std::uint32_t port_count = 0;
};

struct Endpoint {
    bool authored = false;
    std::size_t host = 0;          // generated host index
    std::size_t port = 0;          // region port index (authored endpoint)
    std::size_t cavern = 0;        // region cavern index
    std::string cavern_key, floor_id, label;
    Vec3 position, exit, facing;
    double width = 0, height = 0;
};

struct Built {
    Tunnel tunnel;
    std::vector<Port> new_ports;
    double radius = 0;
};

std::string trim(const std::string& s) {
    std::size_t b = 0, e = s.size();
    while (b < e && std::isspace(static_cast<unsigned char>(s[b]))) ++b;
    while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1]))) --e;
    return s.substr(b, e - b);
}

// Digest of every authored record and reserved box. Generated identities include
// it, so feeding previously generated records back as locked authored content
// can never collide with the identities of the regenerated remainder.
std::string content_digest(const Request& q) {
    Mixer128 m("fbs.caves.authored.v1");
    auto v3 = [&](Vec3 v) { m.f64(v.x).f64(v.y).f64(v.z); };
    auto poly = [&](const std::vector<Vec2>& p) {
        m.u64(p.size());
        for (const auto& v : p) m.f64(v.x).f64(v.y);
    };
    auto strs = [&](const std::vector<std::string>& t) {
        m.u64(t.size());
        for (const auto& x : t) m.str(x);
    };
    auto gen = [&](const Generation& g) {
        m.i64(g.seed).u64(g.parameters.size());
        for (const auto& [k, v] : g.parameters) m.str(k).str(v);
    };
    m.u64(q.reserved.size());
    for (const Bounds& b : q.reserved) { v3(b.minimum); v3(b.maximum); }
    m.u64(q.authored_caverns.size());
    for (const Cavern& c : q.authored_caverns) {
        m.str(c.id).str(c.code).str(c.name);
        v3(c.center); v3(c.bounds.minimum); v3(c.bounds.maximum);
        poly(c.boundary);
        m.u64(static_cast<std::uint64_t>(c.authoring)).u64(c.locked).f64(c.reserved_clearance);
        gen(c.generation);
        strs(c.tags);
    }
    m.u64(q.authored_floors.size());
    for (const Floor& f : q.authored_floors) {
        m.str(f.id).str(f.cavern_id).str(f.code).str(f.name);
        poly(f.boundary);
        m.f64(f.base_z).f64(f.slope_x).f64(f.slope_y).f64(f.variation_amplitude).i64(f.variation_seed).str(f.material);
        m.u64(static_cast<std::uint64_t>(f.traversal));
        strs(f.tags);
    }
    m.u64(q.authored_cliffs.size());
    for (const Cliff& k : q.authored_cliffs) {
        m.str(k.id).str(k.cavern_id).str(k.from_floor_id).str(k.to_floor_id).str(k.code);
        poly(k.boundary);
        m.u64(static_cast<std::uint64_t>(k.transition)).f64(k.height);
    }
    m.u64(q.authored_ports.size());
    for (const Port& p : q.authored_ports) {
        m.str(p.id).str(p.cavern_id).str(p.floor_id).str(p.code).str(p.name);
        v3(p.position); v3(p.facing);
        m.f64(p.width).f64(p.height).u64(static_cast<std::uint64_t>(p.state)).u64(static_cast<std::uint64_t>(p.traversal));
    }
    m.u64(q.authored_tunnels.size());
    for (const Tunnel& t : q.authored_tunnels) {
        m.str(t.id).str(t.start_port_id).str(t.end_port_id).str(t.code).str(t.name).u64(t.centerline.size());
        for (const Profile& p : t.centerline) { v3(p.center); m.f64(p.width).f64(p.height); }
        m.f64(t.maximum_slope_degrees).f64(t.minimum_clearance).u64(static_cast<std::uint64_t>(t.traversal)).u64(t.locked);
        gen(t.generation);
        strs(t.tags);
        m.u64(t.source_node_ids.size());
        for (auto n : t.source_node_ids) m.u64(n);
    }
    const Hash128 h = m.finish();
    char out[40];
    std::snprintf(out, sizeof out, "%016llx%016llx", static_cast<unsigned long long>(h.hi), static_cast<unsigned long long>(h.lo));
    return out;
}

class Generator {
public:
    Generator(const Request& q, const Options& o) : q_(q), o_(o), budget_(o) {}
    Region run();

private:
    // setup
    void setup();
    void copy_authored();
    void add_exclusion(Exclusion e);
    bool clear(Vec3 a, Vec3 b, double threshold, const std::string& exempt);
    bool box_clear(const Bounds& box, double threshold);
    void place_gateways();
    void place_chambers();
    std::size_t add_host_records(Host& h, const Polygon& poly, const Bounds& bounds, Vec3 center, const std::string& role,
                                 std::int32_t seed, Rng& rng);
    void build_lattice();
    // connections
    std::optional<Endpoint> plan_host_port(std::size_t host, Vec2 target, double width, double height, double r);
    std::optional<Endpoint> plan_authored_port(std::size_t port, double r, std::string& why);
    std::optional<Built> route(const Endpoint& a, const Endpoint& b, double r, const std::string& role);
    bool connect_hosts(std::size_t a, std::size_t b, const std::string& role);
    bool connect_authored(std::size_t port, std::size_t host, std::string& why);
    void commit(Built& built, const Endpoint& a, const Endpoint& b);
    void connect_network();
    void connect_authored_ports();
    void finish();
    // Map-scope codes (caverns, tunnels) must stay unique next to authored codes;
    // a colliding generated code gets a deterministic "~n" suffix.
    std::string unique_code(const std::string& base) {
        std::string code = base;
        for (int n = 2; map_codes_.count(fold(code)); ++n) code = base + "~" + std::to_string(n);
        map_codes_[fold(code)] = 1;
        return code;
    }
    static std::string fold(std::string s) {
        for (char& c : s)
            if (c >= 'A' && c <= 'Z') c = static_cast<char>(c - 'A' + 'a');
        return s;
    }
    std::size_t find(std::size_t x) {
        while (parent_[x] != x) x = parent_[x] = parent_[parent_[x]];
        return x;
    }

    const Request& q_;
    const Options& o_;
    Budget budget_;
    Shared sh_;
    Derived dv_;
    std::string sid_, rid_, rtag_;
    std::string content_;  // identity prefix of generated records (region + authored-content digest)
    std::unordered_map<std::string, int> map_codes_;  // map-scope codes already taken
    Region out_;
    std::vector<Exclusion> ex_;
    Grid grid_;
    std::vector<Host> hosts_;
    std::vector<Vec3> lattice_;
    int lx_ = 0, ly_ = 0, lz_ = 0;
    Vec3 lattice_origin_;
    std::vector<std::uint32_t> adj_offset_, adj_;
    std::vector<std::size_t> parent_;
    std::unordered_map<std::string, std::size_t> cavern_index_, port_index_, floor_index_;
    std::vector<bool> port_used_;
    std::uint32_t tunnel_count_ = 0;
    std::uint64_t route_attempts_ = 0, max_route_attempts_ = 0;
    std::vector<std::pair<std::size_t, std::size_t>> host_links_;
    std::map<std::size_t, std::string> failure_;  // last connection failure per host
};

void Generator::setup() {
    sh_ = shared_from_request(q_);
    dv_ = derive(sh_);
    sid_ = settings_id(sh_);
    rid_ = region_identity(sh_, q_.key);
    rtag_ = "R" + key_label(q_.key);
    content_ = rid_ + "/content/" + content_digest(q_);
    out_.version = schema_version;
    out_.algorithm = algorithm_version;
    out_.id = stable_id(rid_);
    out_.settings_id = sid_;
    out_.world_seed = q_.world_seed;
    out_.key = q_.key;
    out_.size = q_.size;
    out_.origin = {static_cast<double>(q_.key.x) * q_.size.x, static_cast<double>(q_.key.y) * q_.size.y,
                   static_cast<double>(q_.key.z) * q_.size.z};
    if (!finite(out_.origin)) throw Error(ErrorCode::invalid_input, "region origin overflows");
    const double cell = std::max({4.0 * dv_.radius, 2.0 * dv_.spacing, std::max(q_.size.x, q_.size.y) / 256.0});
    grid_.init({{-cell, -cell, 0}, {q_.size.x + cell, q_.size.y + cell, q_.size.z}}, cell);
}

void Generator::copy_authored() {
    out_.caverns = q_.authored_caverns;
    out_.floors = q_.authored_floors;
    out_.cliffs = q_.authored_cliffs;
    out_.ports = q_.authored_ports;
    out_.tunnels = q_.authored_tunnels;
    for (std::size_t i = 0; i < out_.caverns.size(); ++i) cavern_index_[uuid_key(out_.caverns[i].id)] = i;
    for (const Cavern& c : out_.caverns) map_codes_[fold(trim(c.code))] = 1;
    for (const Tunnel& t : out_.tunnels) map_codes_[fold(trim(t.code))] = 1;
    for (std::size_t i = 0; i < out_.ports.size(); ++i) port_index_[uuid_key(out_.ports[i].id)] = i;
    for (std::size_t i = 0; i < out_.floors.size(); ++i) floor_index_[uuid_key(out_.floors[i].id)] = i;
    port_used_.assign(out_.ports.size(), false);
    parent_.resize(out_.caverns.size());
    for (std::size_t i = 0; i < parent_.size(); ++i) parent_[i] = i;

    for (const Cavern& c : q_.authored_caverns) {
        Exclusion e;
        e.kind = Exclusion::box;
        e.clearance = c.reserved_clearance;
        e.bounds = c.bounds;
        e.owner = uuid_key(c.id);
        add_exclusion(e);
    }
    for (const Bounds& b : q_.reserved) {
        Exclusion e;
        e.kind = Exclusion::box;
        e.bounds = b;
        add_exclusion(e);
    }
    for (const Tunnel& t : q_.authored_tunnels) {
        const std::size_t a = port_index_.at(uuid_key(t.start_port_id)), b = port_index_.at(uuid_key(t.end_port_id));
        port_used_[a] = port_used_[b] = true;
        const std::size_t ca = cavern_index_.at(uuid_key(out_.ports[a].cavern_id));
        const std::size_t cb = cavern_index_.at(uuid_key(out_.ports[b].cavern_id));
        parent_[find(ca)] = find(cb);
        for (std::size_t k = 1; k < t.centerline.size(); ++k) {
            const Profile& p0 = t.centerline[k - 1];
            const Profile& p1 = t.centerline[k];
            Exclusion e;
            e.kind = Exclusion::capsule;
            e.a = p0.center;
            e.b = p1.center;
            e.radius = section_radius(std::max(p0.width, p1.width), std::max(p0.height, p1.height));
            add_exclusion(e);
        }
    }
}

void Generator::add_exclusion(Exclusion e) {
    if (e.kind == Exclusion::box) {
        e.bounds = inflate(e.bounds, e.clearance);
    } else {
        e.bounds = inflate(segment_bounds(e.a, e.b), e.radius + e.clearance);
    }
    ex_.push_back(e);
    grid_.add(static_cast<std::uint32_t>(ex_.size() - 1), e.bounds);
}

bool Generator::clear(Vec3 a, Vec3 b, double threshold, const std::string& exempt) {
    const Bounds q = inflate(segment_bounds(a, b), threshold);
    return grid_.all(q, [&](std::uint32_t i) {
        const Exclusion& e = ex_[i];
        if (!exempt.empty() && e.owner == exempt) return true;
        if (!boxes_overlap(q, e.bounds)) return true;
        budget_.charge(1);
        return exclusion_distance(e, a, b) >= threshold;
    });
}

bool Generator::box_clear(const Bounds& box, double threshold) {
    const Bounds q = inflate(box, threshold);
    return grid_.all(q, [&](std::uint32_t i) {
        const Exclusion& e = ex_[i];
        if (!boxes_overlap(q, e.bounds)) return true;
        budget_.charge(1);
        const double d = e.kind == Exclusion::box ? box_box_distance(box, inflate(e.bounds, -e.clearance)) - e.clearance
                                                  : segment_box_distance(e.a, e.b, box) - e.radius - e.clearance;
        return d >= threshold;
    });
}

std::size_t Generator::add_host_records(Host& h, const Polygon& poly, const Bounds& bounds, Vec3 center,
                                        const std::string& role, std::int32_t seed, Rng& rng) {
    Cavern c;
    c.id = stable_id(h.identity);
    c.code = unique_code(rtag_ + "-" + h.label);
    c.name = h.name;
    c.center = center;
    c.bounds = bounds;
    c.boundary = poly;
    c.authoring = Authoring::generated;
    c.locked = false;
    c.reserved_clearance = dv_.gap;
    c.generation.seed = seed;
    for (const auto& [k, v] : shared_parameters(sh_)) c.generation.parameters[k] = v;
    c.generation.parameters[generator_key] = algorithm_version;
    c.generation.parameters["fbs.role"] = role;
    c.generation.parameters["fbs.region_id"] = out_.id;
    c.tags = {"generated", role};
    h.cavern = out_.caverns.size();
    h.key = uuid_key(c.id);
    out_.caverns.push_back(c);
    cavern_index_[h.key] = h.cavern;
    parent_.push_back(h.cavern);

    const bool split = h.shelf;
    const std::size_t n = poly.size();
    std::size_t a = 0, b = 0;
    if (split) {
        for (std::size_t i = 0; i < n; ++i)
            if (h.star.angle[i] == h.cut_angle) a = i;
        b = a + n / 2;
    }
    auto add_floor = [&](const Polygon& boundary, double z, const std::string& code, const std::string& name,
                         const std::vector<std::string>& tags) {
        Floor f;
        f.id = stable_id(h.identity + "/floor/" + code);
        f.cavern_id = c.id;
        f.code = code;
        f.name = name;
        f.boundary = boundary;
        f.base_z = z;
        f.variation_seed = static_cast<std::int32_t>(static_cast<std::uint32_t>(rng.next() >> 33));
        f.material = "limestone";
        f.traversal = Traversal::walkable;
        f.tags = tags;
        floor_index_[uuid_key(f.id)] = out_.floors.size();
        out_.floors.push_back(f);
        return out_.floors.size() - 1;
    };
    if (!split) {
        h.main_floor = add_floor(poly, h.floor_z, "F1", h.name + " floor", {"generated"});
        return h.cavern;
    }
    Polygon main, shelf;
    for (std::size_t i = a; i <= b; ++i) main.push_back(poly[i]);
    for (std::size_t i = b; i < n; ++i) shelf.push_back(poly[i]);
    for (std::size_t i = 0; i <= a; ++i) shelf.push_back(poly[i]);
    h.main_floor = add_floor(main, h.floor_z, "F1", h.name + " floor", {"generated"});
    h.shelf_floor = add_floor(shelf, h.floor_z + h.shelf_height, "F2", h.name + " raised shelf", {"generated", "shelf"});
    Cliff k;
    k.id = stable_id(h.identity + "/cliff/K1");
    k.cavern_id = c.id;
    k.from_floor_id = out_.floors[h.main_floor].id;
    k.to_floor_id = out_.floors[h.shelf_floor].id;
    k.code = "K1";
    k.boundary = {poly[a], poly[b]};
    k.transition = h.shelf_height > 1.0 ? Transition::climbable_cliff : Transition::stairs;
    k.height = h.shelf_height;
    out_.cliffs.push_back(k);
    return h.cavern;
}

void Generator::place_gateways() {
    for (int f = 0; f < 6; ++f) {
        if (!(q_.face_mask & (1u << f))) continue;
        budget_.check_cancelled();
        const Face face = static_cast<Face>(f);
        const Seam seam = derive_seam(sh_, dv_, q_.key, face);
        Host h;
        h.gateway = true;
        h.identity = content_ + "/gateway/" + face_label(face);
        h.label = std::string("G-") + face_label(face);
        h.name = std::string("Gateway ") + face_label(face);
        Mixer128 mixer("fbs.caves.gateway.v1");
        mixer.str(rid_ + "/gateway/" + face_label(face)); // shape independent of authored content
        const Hash128 hash = mixer.finish();
        Rng rng(hash);
        const double psi = std::atan2(seam.facing.y, seam.facing.x);
        h.star = make_star(rng, {0, 0}, dv_.gateway_radius, 12);
        const double rb = h.star.boundary(psi);
        const Vec2 center{seam.port.x - std::cos(psi) * (rb - port_inset), seam.port.y - std::sin(psi) * (rb - port_inset)};
        h.star.center = center;
        h.floor_z = seam.port.z - 0.5 * seam.height;
        h.ceiling_z = h.floor_z + dv_.gateway_height;
        const Polygon poly = h.star.polygon();
        const Bounds bounds = polygon_bounds(poly, h.floor_z, h.ceiling_z);
        const double m = 0.5 * dv_.gap;
        if (bounds.minimum.x < m || bounds.minimum.y < m || bounds.minimum.z < m || bounds.maximum.x > q_.size.x - m ||
            bounds.maximum.y > q_.size.y - m || bounds.maximum.z > q_.size.z - m)
            throw Error(ErrorCode::constraint, std::string("gateway for face ") + face_label(face) + " does not fit inside the region");
        if (!box_clear(bounds, 2.0 * dv_.radius + dv_.gap))
            throw Error(ErrorCode::constraint, std::string("gateway for face ") + face_label(face) +
                                                   " conflicts with authored, reserved or other gateway content");
        const double rs = section_radius(seam.width, seam.height);
        if (!clear(seam.port, seam.position, rs + dv_.gap, ""))
            throw Error(ErrorCode::constraint, std::string("seam corridor for face ") + face_label(face) +
                                                   " is blocked by authored, reserved or gateway content");
        const Vec3 c3{center.x, center.y, 0.5 * (h.floor_z + h.ceiling_z)};
        add_host_records(h, poly, bounds, c3, "gateway", seed32(hash), rng);
        Exclusion box;
        box.kind = Exclusion::box;
        box.bounds = bounds;
        box.owner = h.key;
        add_exclusion(box);
        Exclusion corridor;
        corridor.kind = Exclusion::capsule;
        corridor.a = seam.port;
        corridor.b = seam.position;
        corridor.radius = rs;
        corridor.owner = "corridor:" + h.key;
        add_exclusion(corridor);

        Port p;
        p.id = stable_id(h.identity + "/port/seam");
        p.cavern_id = out_.caverns[h.cavern].id;
        p.floor_id = out_.floors[h.main_floor].id;
        p.code = "P0";
        p.name = h.name + " seam";
        p.position = seam.port;
        p.facing = seam.facing;
        p.width = seam.width;
        p.height = seam.height;
        p.state = PortState::optional;
        p.traversal = seam.slope_degrees <= walkable_slope_degrees ? Traversal::walkable : Traversal::difficult;
        if (!polygon_contains(out_.floors[h.main_floor].boundary, {p.position.x, p.position.y}))
            throw Error(ErrorCode::internal, "seam port is outside its gateway floor");
        port_index_[uuid_key(p.id)] = out_.ports.size();
        out_.ports.push_back(p);
        port_used_.push_back(false);
        h.ports.push_back(seam.port);
        h.port_radius.push_back(rs);

        Boundary bd;
        bd.key = seam.key;
        bd.port_id = p.id;
        bd.face = face;
        bd.position = seam.position;
        bd.width = seam.width;
        bd.height = seam.height;
        out_.boundaries.push_back(bd);
        hosts_.push_back(std::move(h));
    }
}

void Generator::place_chambers() {
    if (q_.chamber_count == 0) return;
    const double zlo = dv_.radius + dv_.gap, zhi = q_.size.z - q_.chamber_height - dv_.gap;
    if (zhi < zlo) throw Error(ErrorCode::constraint, "chamber_height does not fit inside the region height");
    if (q_.tunnel_height + 0.5 > q_.chamber_height)
        throw Error(ErrorCode::constraint, "chamber_height must exceed tunnel_height by at least 0.5 m");
    // Each passage port needs 2*radius + gap of chamber wall; the smallest
    // chamber must be able to host at least three (an in/out pair plus a branch).
    const double port_spacing = 2.0 * dv_.radius + dv_.gap;
    const double capacity = std::floor(2.0 * pi * 0.72 * q_.chamber_radius_min / port_spacing);
    if (capacity < 3)
        throw Error(ErrorCode::constraint, "chamber_radius_min " + format_double(q_.chamber_radius_min) +
                                               " m is too small to host three passage ports of this size (needs about " +
                                               format_double(std::ceil(3.0 * port_spacing / (2.0 * pi * 0.72))) + " m)");
    const double edge = q_.chamber_radius_min + dv_.gap;
    if (2.0 * edge >= std::min(q_.size.x, q_.size.y))
        throw Error(ErrorCode::constraint, "chamber_radius_min does not fit inside the region");
    Mixer128 mixer("fbs.caves.chambers.v1");
    mixer.str(rid_);
    Rng rng(mixer.finish());
    const std::uint64_t attempts = 200ull * q_.chamber_count;
    std::uint32_t placed = 0;
    for (std::uint64_t attempt = 0; attempt < attempts && placed < q_.chamber_count; ++attempt) {
        budget_.charge(8);
        const double r = rng.range(q_.chamber_radius_min, q_.chamber_radius_max);
        const std::size_t n = 2 * (5 + rng.below(4));
        const double mx = r + dv_.gap;
        if (2.0 * mx >= std::min(q_.size.x, q_.size.y)) continue;
        const Vec2 center{rng.range(mx, q_.size.x - mx), rng.range(mx, q_.size.y - mx)};
        const double floor_z = rng.range(zlo, zhi);
        Rng shape(Hash128{rng.next(), rng.next()});
        Star star = make_star(shape, center, r, n);
        const Polygon poly = star.polygon();
        const Bounds bounds = polygon_bounds(poly, floor_z, floor_z + q_.chamber_height);
        // Room for two facing access stubs (each ends r + gap outside its box) plus routing space.
        if (!box_clear(bounds, 2.0 * (dv_.radius + dv_.gap) + std::max(2.0, 0.5 * dv_.spacing))) continue;

        Host h;
        h.identity = content_ + "/chamber/" + std::to_string(placed);
        char label[16];
        std::snprintf(label, sizeof label, "C%02u", placed + 1);
        h.label = label;
        h.name = "Chamber " + std::to_string(placed + 1);
        h.star = star;
        h.floor_z = floor_z;
        h.ceiling_z = floor_z + q_.chamber_height;
        const double shelf_hi = std::min(3.0, q_.chamber_height - dv_.height_max - 0.75);
        const bool roomy = r >= 1.5 * q_.tunnel_width + 2.0;
        if (q_.shelves && shape.uniform() < 0.5 && shelf_hi >= 0.8 && roomy) {
            h.shelf = true;
            h.shelf_height = shape.range(0.8, shelf_hi);
            h.cut_angle = star.angle[shape.below(static_cast<std::uint32_t>(n / 2))];
        }
        Mixer128 hm("fbs.caves.chamber.v1");
        hm.str(h.identity);
        const Hash128 hash = hm.finish();
        const Vec3 c3{center.x, center.y, floor_z + 0.5 * q_.chamber_height};
        add_host_records(h, poly, bounds, c3, "chamber", seed32(hash), shape);
        Exclusion e;
        e.kind = Exclusion::box;
        e.bounds = bounds;
        e.owner = h.key;
        add_exclusion(e);
        hosts_.push_back(std::move(h));
        ++placed;
    }
    if (placed < q_.chamber_count)
        throw Error(ErrorCode::constraint, "only " + std::to_string(placed) + " of " + std::to_string(q_.chamber_count) +
                                               " chambers fit without overlapping existing or reserved content");
}

void Generator::build_lattice() {
    const double s = dv_.spacing;
    const double m = dv_.radius + 0.5 * dv_.gap;
    auto count = [&](double extent) {
        const double usable = extent - 2.0 * m;
        return usable < 0 ? 0.0 : std::floor(usable / s) + 1.0;
    };
    const double cx = count(q_.size.x), cy = count(q_.size.y), cz = count(q_.size.z);
    const double total = cx * cy * cz;
    const double cap = static_cast<double>(std::min<std::uint64_t>(o_.max_points, hard_max_support_points));
    if (total > cap)
        throw Error(ErrorCode::budget, "support lattice needs " + format_double(total) + " points; limit is " + format_double(cap));
    lx_ = static_cast<int>(cx);
    ly_ = static_cast<int>(cy);
    lz_ = static_cast<int>(cz);
    if (total == 0) return;
    Mixer128 mixer("fbs.caves.lattice.v1");
    mixer.str(rid_);
    Rng rng(mixer.finish());
    auto offset = [&](double extent, int n) { return m + 0.5 * ((extent - 2.0 * m) - (n - 1) * s); };
    const double ox = offset(q_.size.x, lx_), oy = offset(q_.size.y, ly_), oz = offset(q_.size.z, lz_);
    lattice_origin_ = {ox, oy, oz};
    lattice_.reserve(static_cast<std::size_t>(total));
    for (int k = 0; k < lz_; ++k)
        for (int j = 0; j < ly_; ++j)
            for (int i = 0; i < lx_; ++i) {
                const double jx = rng.range(-0.3, 0.3) * s, jy = rng.range(-0.3, 0.3) * s, jz = rng.range(-0.3, 0.3) * s;
                lattice_.push_back({std::clamp(ox + i * s + jx, m, q_.size.x - m), std::clamp(oy + j * s + jy, m, q_.size.y - m),
                                    std::clamp(oz + k * s + jz, m, q_.size.z - m)});
            }
    budget_.charge(lattice_.size());

    const std::uint64_t edge_cap = std::min<std::uint64_t>(o_.max_edges, hard_max_support_edges);
    const std::uint64_t k = q_.neighbors;
    if (static_cast<double>(lattice_.size()) * static_cast<double>(k) > static_cast<double>(edge_cap))
        throw Error(ErrorCode::budget, "support graph would exceed the edge limit");
    const double radius = 2.2 * s;
    const double slope_limit = q_.maximum_slope_degrees - slope_margin_degrees;
    std::vector<std::pair<std::uint32_t, std::uint32_t>> pairs;
    pairs.reserve(lattice_.size() * k);
    std::vector<std::pair<double, std::uint32_t>> cand;
    auto idx = [&](int i, int j, int kk) { return static_cast<std::uint32_t>((kk * ly_ + j) * lx_ + i); };
    for (int kk = 0; kk < lz_; ++kk)
        for (int j = 0; j < ly_; ++j)
            for (int i = 0; i < lx_; ++i) {
                const std::uint32_t u = idx(i, j, kk);
                cand.clear();
                for (int dk = -2; dk <= 2; ++dk)
                    for (int dj = -2; dj <= 2; ++dj)
                        for (int di = -2; di <= 2; ++di) {
                            const int a = i + di, b = j + dj, c = kk + dk;
                            if ((di | dj | dk) == 0 || a < 0 || b < 0 || c < 0 || a >= lx_ || b >= ly_ || c >= lz_) continue;
                            const std::uint32_t v = idx(a, b, c);
                            const double d = distance(lattice_[u], lattice_[v]);
                            if (d > radius) continue;
                            if (segment_slope_degrees(lattice_[u], lattice_[v]) > slope_limit) continue;
                            cand.push_back({d, v});
                        }
                budget_.charge(125);
                std::sort(cand.begin(), cand.end());
                for (std::size_t c = 0; c < cand.size() && c < k; ++c)
                    pairs.push_back({std::min(u, cand[c].second), std::max(u, cand[c].second)});
            }
    std::sort(pairs.begin(), pairs.end());
    pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());
    adj_offset_.assign(lattice_.size() + 1, 0);
    for (const auto& [a, b] : pairs) {
        ++adj_offset_[a + 1];
        ++adj_offset_[b + 1];
    }
    for (std::size_t i = 1; i < adj_offset_.size(); ++i) adj_offset_[i] += adj_offset_[i - 1];
    adj_.assign(adj_offset_.back(), 0);
    std::vector<std::uint32_t> fill(adj_offset_.begin(), adj_offset_.end() - 1);
    for (const auto& [a, b] : pairs) {
        adj_[fill[a]++] = b;
        adj_[fill[b]++] = a;
    }
    budget_.charge(pairs.size());
}

std::optional<Endpoint> Generator::plan_host_port(std::size_t hi, Vec2 target, double width, double height, double r) {
    Host& h = hosts_[hi];
    const double psi0 = std::atan2(target.y - h.star.center.y, target.x - h.star.center.x);
    const double step = 7.5 * pi / 180.0;
    const double m = r + 0.5 * dv_.gap;
    const Bounds& cb = out_.caverns[h.cavern].bounds;
    for (int attempt = 0; attempt < 48; ++attempt) {
        const int k = (attempt + 1) / 2;
        const double psi = psi0 + (attempt % 2 == 1 ? k : -k) * step;
        const double rho = h.star.boundary(psi) - port_inset;
        if (rho <= 0.5) continue;
        const Vec2 dir{std::cos(psi), std::sin(psi)};
        const Vec2 xy{h.star.center.x + dir.x * rho, h.star.center.y + dir.y * rho};
        bool shelf = false;
        if (h.shelf) {
            const double rel = norm_angle(psi - h.cut_angle);
            shelf = rel > pi;
            if (rho * std::fabs(std::sin(rel)) < 0.5 * width + 0.5) continue;
        }
        const double floor_z = h.floor_z + (shelf ? h.shelf_height : 0.0);
        const Vec3 pos{xy.x, xy.y, floor_z + 0.5 * height};
        if (pos.z + 0.5 * height > h.ceiling_z - 0.25) continue;
        const std::size_t fl = shelf ? h.shelf_floor : h.main_floor;
        if (!polygon_contains(out_.floors[fl].boundary, xy)) continue;
        bool separated = true;
        for (std::size_t i = 0; i < h.ports.size() && separated; ++i)
            separated = distance(h.ports[i], pos) >= h.port_radius[i] + r + dv_.gap;
        if (!separated) continue;
        // Horizontal access stub to just outside the chamber box inflated by r + gap.
        const double inf = r + dv_.gap;
        double t = std::numeric_limits<double>::infinity();
        if (dir.x > 1e-12) t = std::min(t, (cb.maximum.x + inf - xy.x) / dir.x);
        if (dir.x < -1e-12) t = std::min(t, (cb.minimum.x - inf - xy.x) / dir.x);
        if (dir.y > 1e-12) t = std::min(t, (cb.maximum.y + inf - xy.y) / dir.y);
        if (dir.y < -1e-12) t = std::min(t, (cb.minimum.y - inf - xy.y) / dir.y);
        if (!std::isfinite(t) || t <= 0) continue;
        const Vec3 exit{xy.x + dir.x * (t + 0.05), xy.y + dir.y * (t + 0.05), pos.z};
        if (exit.x < m || exit.y < m || exit.z < m || exit.x > q_.size.x - m || exit.y > q_.size.y - m || exit.z > q_.size.z - m)
            continue;
        if (!clear(pos, exit, r + dv_.gap, h.key)) continue;
        Endpoint e;
        e.host = hi;
        e.cavern = h.cavern;
        e.cavern_key = h.key;
        e.floor_id = out_.floors[fl].id;
        e.position = pos;
        e.exit = exit;
        if (!stable_unit({dir.x, dir.y, 0.0}, e.facing)) continue;
        e.width = width;
        e.height = height;
        e.label = h.label;
        return e;
    }
    return std::nullopt;
}

std::optional<Endpoint> Generator::plan_authored_port(std::size_t pi_, double r, std::string& why) {
    const Port& p = out_.ports[pi_];
    const std::size_t cv = cavern_index_.at(uuid_key(p.cavern_id));
    const Cavern& c = out_.caverns[cv];
    const double hl = std::hypot(p.facing.x, p.facing.y);
    if (hl < 1e-6) {
        why = "port " + p.id + " faces vertically; generated passages leave ports horizontally";
        return std::nullopt;
    }
    const Vec2 dir{p.facing.x / hl, p.facing.y / hl};
    const double inf = r + dv_.gap + c.reserved_clearance;
    double t = std::numeric_limits<double>::infinity();
    if (dir.x > 1e-12) t = std::min(t, (c.bounds.maximum.x + inf - p.position.x) / dir.x);
    if (dir.x < -1e-12) t = std::min(t, (c.bounds.minimum.x - inf - p.position.x) / dir.x);
    if (dir.y > 1e-12) t = std::min(t, (c.bounds.maximum.y + inf - p.position.y) / dir.y);
    if (dir.y < -1e-12) t = std::min(t, (c.bounds.minimum.y - inf - p.position.y) / dir.y);
    const Vec3 exit{p.position.x + dir.x * (t + 0.05), p.position.y + dir.y * (t + 0.05), p.position.z};
    const double m = r + 0.5 * dv_.gap;
    if (!std::isfinite(t) || exit.x < m || exit.y < m || exit.z < m || exit.x > q_.size.x - m || exit.y > q_.size.y - m ||
        exit.z > q_.size.z - m) {
        why = "port " + p.id + " access stub would leave the region";
        return std::nullopt;
    }
    if (!clear(p.position, exit, r + dv_.gap, uuid_key(c.id))) {
        why = "port " + p.id + " access stub is blocked by other content";
        return std::nullopt;
    }
    Endpoint e;
    e.authored = true;
    e.port = pi_;
    e.cavern = cv;
    e.cavern_key = uuid_key(c.id);
    e.floor_id = p.floor_id;
    e.position = p.position;
    e.exit = exit;
    e.facing = {dir.x, dir.y, 0.0};
    e.width = p.width;
    e.height = p.height;
    e.label = p.code;
    return e;
}

std::optional<Built> Generator::route(const Endpoint& A, const Endpoint& B, double r, const std::string& role) {
    if (lattice_.empty()) return std::nullopt;
    const double s = dv_.spacing;
    const double threshold = r + dv_.gap;
    const double m = r + 0.5 * dv_.gap;
    const double slope_limit = q_.maximum_slope_degrees - slope_margin_degrees;
    const double L = distance(A.exit, B.exit);
    if (L < keypoint_spacing) return std::nullopt;
    const std::uint32_t tunnel_index = tunnel_count_;

    for (int attempt = 0; attempt < 2; ++attempt) {
        if (route_attempts_ >= max_route_attempts_) return std::nullopt;
        ++route_attempts_;
        budget_.check_cancelled();
        const double sum = attempt == 0 ? L * 1.35 + 4.0 * s : L * 2.2 + 10.0 * s;
        const Vec3 mid = mul(add(A.exit, B.exit), 0.5);
        const double half = 0.5 * sum;
        // Lattice index window covering the ellipse (node i lies within
        // origin + i*s +- 0.3*s, then clamped into the region margin).
        auto range = [&](double c, double origin, int n) {
            const double lo_f = std::floor((c - half - origin) / s) - 1.0;
            const double hi_f = std::ceil((c + half - origin) / s) + 1.0;
            const int lo = lo_f < 0 ? 0 : (lo_f > n - 1 ? n - 1 : static_cast<int>(lo_f));
            const int hi = hi_f < 0 ? -1 : (hi_f > n - 1 ? n - 1 : static_cast<int>(hi_f));
            return std::pair<int, int>(lo, hi);
        };
        const auto [i0, i1] = range(mid.x, lattice_origin_.x, lx_);
        const auto [j0, j1] = range(mid.y, lattice_origin_.y, ly_);
        const auto [k0, k1] = range(mid.z, lattice_origin_.z, lz_);

        RouteGraph g;
        g.nodes.reserve(64);
        const bool a_high = A.exit.z >= B.exit.z;
        g.nodes.push_back(a_high ? A.exit : B.exit);
        g.nodes.push_back(a_high ? B.exit : A.exit);
        std::unordered_map<std::uint32_t, std::uint32_t> local;
        std::vector<std::uint32_t> global;
        const std::uint64_t point_cap = std::min<std::uint64_t>(o_.max_points, hard_max_support_points);
        for (int k = k0; k <= k1; ++k)
            for (int j = j0; j <= j1; ++j)
                for (int i = i0; i <= i1; ++i) {
                    const std::uint32_t id = static_cast<std::uint32_t>((k * ly_ + j) * lx_ + i);
                    const Vec3 p = lattice_[id];
                    budget_.charge(1);
                    if (distance(p, A.exit) + distance(p, B.exit) > sum) continue;
                    if (p.x < m || p.y < m || p.z < m || p.x > q_.size.x - m || p.y > q_.size.y - m || p.z > q_.size.z - m)
                        continue;
                    if (distance(p, A.exit) < keypoint_spacing || distance(p, B.exit) < keypoint_spacing) continue;
                    if (!clear(p, p, threshold, "")) continue;
                    if (g.nodes.size() + 1 > point_cap)
                        throw Error(ErrorCode::budget, "route support subgraph exceeds max_points");
                    local.emplace(id, static_cast<std::uint32_t>(g.nodes.size()));
                    global.push_back(id);
                    g.nodes.push_back(p);
                }
        for (std::size_t li = 0; li < global.size(); ++li) {
            const std::uint32_t u = global[li];
            for (std::uint32_t e = adj_offset_[u]; e < adj_offset_[u + 1]; ++e) {
                const std::uint32_t v = adj_[e];
                if (v <= u) continue;
                const auto it = local.find(v);
                if (it == local.end()) continue;
                if (!clear(lattice_[u], lattice_[v], threshold, "")) continue;
                g.edges.push_back({static_cast<std::uint32_t>(li + 2), it->second});
            }
            if (g.edges.size() > std::min<std::uint64_t>(o_.max_edges, hard_max_support_edges))
                throw Error(ErrorCode::budget, "route support subgraph exceeds max_edges");
        }
        // Access points join their nearest admissible support nodes (and each other).
        for (std::uint32_t key = 0; key < 2; ++key) {
            const Vec3 kp = g.nodes[key];
            std::vector<std::pair<double, std::uint32_t>> near;
            for (std::size_t li = 0; li < global.size(); ++li) {
                const double d = distance(kp, g.nodes[li + 2]);
                if (d <= 2.2 * s) near.push_back({d, static_cast<std::uint32_t>(li + 2)});
            }
            std::sort(near.begin(), near.end());
            std::uint32_t taken = 0;
            for (const auto& [d, v] : near) {
                if (taken >= q_.neighbors) break;
                if (segment_slope_degrees(kp, g.nodes[v]) > slope_limit) continue;
                if (!clear(kp, g.nodes[v], threshold, "")) continue;
                g.edges.push_back({key, v});
                ++taken;
            }
        }
        if (L <= 2.2 * s && segment_slope_degrees(A.exit, B.exit) <= slope_limit && clear(A.exit, B.exit, threshold, ""))
            g.edges.push_back({0, 1});

        // Cheap connectivity pre-check before invoking the engine.
        {
            std::vector<std::vector<std::uint32_t>> adjacency(g.nodes.size());
            for (const auto& [a, b] : g.edges) {
                adjacency[a].push_back(b);
                adjacency[b].push_back(a);
            }
            std::vector<bool> seen(g.nodes.size(), false);
            std::vector<std::uint32_t> stack{0};
            seen[0] = true;
            while (!stack.empty()) {
                const std::uint32_t u = stack.back();
                stack.pop_back();
                for (std::uint32_t v : adjacency[u])
                    if (!seen[v]) {
                        seen[v] = true;
                        stack.push_back(v);
                    }
            }
            budget_.charge(g.nodes.size() + g.edges.size());
            if (!seen[1]) continue;
        }

        g.domain = {{0, 0, 0}, q_.size};
        Mixer128 rm("fbs.caves.route.v1");
        rm.str(rid_).u64(tunnel_index).u64(route_attempts_).u64(static_cast<std::uint64_t>(attempt));
        const Hash128 rh = rm.finish();
        RouteSettings settings;
        settings.seed = std::max<std::int32_t>(1, seed32(rh));
        settings.neighbors = q_.neighbors;
        // The engine runs on what is left of the one regional work budget and its
        // actual counter is charged back, including for failed or aborted runs.
        if (budget_.remaining() == 0) throw Error(ErrorCode::budget, "regional work limit exceeded");
        std::uint64_t native_work = 0;
        RouteLimits limits;
        limits.work_limit = budget_.remaining();
        limits.max_points = o_.max_points;
        limits.max_edges = o_.max_edges;
        limits.cancelled = &o_.cancelled;
        limits.work_used = &native_work;
        const RouteResult rr = route_karst(g, settings, limits);
        budget_.charge(native_work);
        if (!rr.found) continue;

        std::vector<std::uint32_t> path = rr.path;
        std::vector<std::uint32_t> native = rr.native_node_ids;
        if (!a_high) {
            std::reverse(path.begin(), path.end());
            std::reverse(native.begin(), native.end());
        }
        // Profiles: exact port endpoints, stub exits with port dimensions,
        // interior widths/heights varied deterministically.
        Rng prng(rh);
        Tunnel t;
        const double h_t = A.height;
        t.centerline.push_back({A.position, A.width, A.height});
        for (std::size_t i = 0; i < path.size(); ++i) {
            const Vec3 c = g.nodes[path[i]];
            if (i == 0) {
                t.centerline.push_back({c, A.width, h_t});
            } else if (i + 1 == path.size()) {
                t.centerline.push_back({c, B.width, h_t});
            } else {
                const double w = q_.tunnel_width * (1.0 + 0.3 * prng.uniform());
                const double hh = h_t == q_.tunnel_height ? q_.tunnel_height * (1.0 + 0.15 * prng.uniform()) : h_t;
                t.centerline.push_back({c, w, hh});
            }
        }
        t.centerline.push_back({B.position, B.width, B.height});
        // Floor slope must respect the limit as well as the centreline slope.
        const double floor_limit = q_.maximum_slope_degrees - 0.5 * slope_margin_degrees;
        for (std::size_t pass = 0; pass < t.centerline.size(); ++pass) {
            bool changed = false;
            for (std::size_t i = 2; i + 2 < t.centerline.size(); ++i) {
                auto floor_slope = [&](std::size_t a, std::size_t b) {
                    const Profile& pa = t.centerline[a];
                    const Profile& pb = t.centerline[b];
                    return segment_slope_degrees({pa.center.x, pa.center.y, pa.center.z - 0.5 * pa.height},
                                                 {pb.center.x, pb.center.y, pb.center.z - 0.5 * pb.height});
                };
                if (t.centerline[i].height != h_t && (floor_slope(i - 1, i) > floor_limit || floor_slope(i, i + 1) > floor_limit)) {
                    t.centerline[i].height = h_t;
                    changed = true;
                }
            }
            if (!changed) break;
        }
        // Final verification of the whole swept passage.
        bool ok = true;
        double max_slope = 0;
        for (std::size_t i = 1; i < t.centerline.size() && ok; ++i) {
            const Vec3 a = t.centerline[i - 1].center, b = t.centerline[i].center;
            const double sl = segment_slope_degrees(a, b);
            max_slope = std::max(max_slope, sl);
            if (sl > slope_limit) ok = false;
            const std::string exempt = i == 1 ? A.cavern_key : (i + 1 == t.centerline.size() ? B.cavern_key : "");
            if (ok && !clear(a, b, threshold, exempt)) ok = false;
        }
        if (!ok) continue;

        Built built;
        built.radius = r;
        t.start_port_id = "";
        t.end_port_id = "";
        t.maximum_slope_degrees = q_.maximum_slope_degrees;
        t.minimum_clearance = h_t;
        for (const Profile& p : t.centerline) t.minimum_clearance = std::min(t.minimum_clearance, p.height);
        t.traversal = max_slope <= walkable_slope_degrees ? Traversal::walkable : Traversal::difficult;
        t.locked = false;
        t.generation.seed = settings.seed;
        for (const auto& [k, v] : shared_parameters(sh_)) t.generation.parameters[k] = v;
        t.generation.parameters[generator_key] = algorithm_version;
        t.generation.parameters["fbs.role"] = role;
        t.generation.parameters["fbs.region_id"] = out_.id;
        t.generation.parameters["fbs.route.engine"] = "KarstNSim::run_simulation_memory";
        t.generation.parameters["fbs.route.seed"] = std::to_string(settings.seed);
        t.generation.parameters["fbs.route.support_nodes"] = std::to_string(g.nodes.size());
        t.generation.parameters["fbs.route.support_edges"] = std::to_string(g.edges.size());
        t.generation.parameters["fbs.route.source_node_ids"] = "native skeleton node ids of this route run";
        t.tags = {"generated", role};
        t.source_node_ids = native;
        built.tunnel = std::move(t);
        return built;
    }
    return std::nullopt;
}

void Generator::commit(Built& built, const Endpoint& A, const Endpoint& B) {
    auto make_port = [&](const Endpoint& e) -> std::string {
        if (e.authored) {
            port_used_[e.port] = true;
            return out_.ports[e.port].id;
        }
        Host& h = hosts_[e.host];
        ++h.port_count;
        Port p;
        p.id = stable_id(h.identity + "/port/" + std::to_string(h.port_count));
        p.cavern_id = out_.caverns[h.cavern].id;
        p.floor_id = e.floor_id;
        p.code = "P" + std::to_string(h.port_count);
        p.name = h.name + " passage " + std::to_string(h.port_count);
        p.position = e.position;
        p.facing = e.facing;
        p.width = e.width;
        p.height = e.height;
        p.state = PortState::generated;
        p.traversal = Traversal::walkable;
        port_index_[uuid_key(p.id)] = out_.ports.size();
        out_.ports.push_back(p);
        port_used_.push_back(true);
        h.ports.push_back(e.position);
        h.port_radius.push_back(built.radius);
        return p.id;
    };
    Tunnel& t = built.tunnel;
    t.start_port_id = make_port(A);
    t.end_port_id = make_port(B);
    ++tunnel_count_;
    t.id = stable_id(content_ + "/tunnel/" + uuid_key(t.start_port_id) + "/" + uuid_key(t.end_port_id));
    char code[32];
    std::snprintf(code, sizeof code, "-T%03u", tunnel_count_);
    t.code = unique_code(rtag_ + code);
    t.name = "Passage " + A.label + " to " + B.label;
    for (std::size_t i = 1; i < t.centerline.size(); ++i) {
        Exclusion e;
        e.kind = Exclusion::capsule;
        e.a = t.centerline[i - 1].center;
        e.b = t.centerline[i].center;
        e.radius = built.radius;
        add_exclusion(e);
    }
    parent_[find(A.cavern)] = find(B.cavern);
    out_.tunnels.push_back(std::move(t));
}

bool Generator::connect_hosts(std::size_t a, std::size_t b, const std::string& role) {
    const double r = dv_.radius;
    const Host& ha = hosts_[a];
    const Host& hb = hosts_[b];
    auto A = plan_host_port(a, hb.star.center, q_.tunnel_width, q_.tunnel_height, r);
    if (!A) {
        failure_[a] = "no room for another admissible port and access stub on " + ha.name;
        return false;
    }
    auto B = plan_host_port(b, ha.star.center, q_.tunnel_width, q_.tunnel_height, r);
    if (!B) {
        failure_[b] = "no room for another admissible port and access stub on " + hb.name;
        return false;
    }
    auto built = route(*A, *B, r, role);
    if (!built) {
        failure_[a] = failure_[b] = "no clear slope-limited route between " + ha.name + " and " + hb.name;
        return false;
    }
    commit(*built, *A, *B);
    host_links_.push_back({a, b});
    return true;
}

bool Generator::connect_authored(std::size_t port, std::size_t host, std::string& why) {
    const Port& p = out_.ports[port];
    const double r = section_radius(std::max(dv_.width_max, p.width), std::max(dv_.height_max, p.height));
    auto A = plan_authored_port(port, r, why);
    if (!A) return false;
    auto B = plan_host_port(host, {A->exit.x, A->exit.y}, p.width, p.height, r);
    if (!B) {
        why = "no admissible port on " + hosts_[host].name + " for the dimensions of port " + p.id;
        return false;
    }
    auto built = route(*A, *B, r, "authored-link");
    if (!built) {
        why = "no admissible route from port " + p.id + " to " + hosts_[host].name;
        return false;
    }
    commit(*built, *A, *B);
    return true;
}

void Generator::connect_network() {
    const std::size_t n = hosts_.size();
    if (n == 0) {
        if (q_.loop_count > 0) throw Error(ErrorCode::constraint, "loops requested but the region has no chambers or gateways");
        return;
    }
    const std::uint64_t possible = n * (n - 1) / 2 - (n - 1);
    if (q_.loop_count > possible)
        throw Error(ErrorCode::constraint, std::to_string(q_.loop_count) + " loops requested but only " +
                                               std::to_string(possible) + " distinct loop connections exist");
    std::vector<std::tuple<double, std::size_t, std::size_t>> pairs;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = i + 1; j < n; ++j) {
            const Vec3 a = out_.caverns[hosts_[i].cavern].center, b = out_.caverns[hosts_[j].cavern].center;
            pairs.emplace_back(distance(a, b), i, j);
        }
    std::sort(pairs.begin(), pairs.end());
    budget_.charge(pairs.size());
    // Spanning network: repeatedly take the cheapest untried pair joining two
    // components, where cost is distance inflated by the hosts' current degree.
    // The degree term spreads passages over chambers so wall capacity for ports
    // is not exhausted by a few hubs (plain shortest-first can strand chambers).
    std::vector<std::size_t> comp(n);
    for (std::size_t i = 0; i < n; ++i) comp[i] = i;
    auto root = [&](std::size_t x) {
        while (comp[x] != x) x = comp[x] = comp[comp[x]];
        return x;
    };
    std::size_t components = n;
    std::vector<bool> used_pair(pairs.size(), false), tried(pairs.size(), false);
    std::vector<std::uint32_t> degree(n, 0);
    while (components > 1 && route_attempts_ < max_route_attempts_) {
        std::size_t best = pairs.size();
        double best_cost = std::numeric_limits<double>::infinity();
        for (std::size_t p = 0; p < pairs.size(); ++p) {
            if (tried[p]) continue;
            const auto& [d, a, b] = pairs[p];
            if (root(a) == root(b)) continue;
            const double cost = d * (1.0 + 0.35 * (degree[a] + degree[b]));
            if (cost < best_cost) {
                best_cost = cost;
                best = p;
            }
        }
        budget_.charge(pairs.size());
        if (best == pairs.size()) break;
        tried[best] = true;
        const auto [d, a, b] = pairs[best];
        if (connect_hosts(a, b, "passage")) {
            comp[root(a)] = root(b);
            --components;
            used_pair[best] = true;
            ++degree[a];
            ++degree[b];
        }
    }
    if (components > 1) {
        for (std::size_t i = 0; i < n; ++i)
            if (root(i) != root(0)) {
                const std::string why = failure_.count(i) ? failure_[i]
                                        : route_attempts_ >= max_route_attempts_ ? "route attempt limit reached"
                                                                                 : "no candidate connection";
                throw Error(ErrorCode::no_route, hosts_[i].name + " could not be connected to the regional network (" + why + ")");
            }
    }
    // Deliberate loops: extra passages between chambers not directly linked.
    std::uint32_t loops = 0;
    for (std::size_t p = 0; p < pairs.size() && loops < q_.loop_count; ++p) {
        if (used_pair[p]) continue;
        const auto [d, a, b] = pairs[p];
        bool direct = false;
        for (const auto& [x, y] : host_links_) direct = direct || (x == a && y == b) || (x == b && y == a);
        if (direct) continue;
        if (route_attempts_ >= max_route_attempts_) break;
        if (connect_hosts(a, b, "loop")) {
            used_pair[p] = true;
            ++loops;
        }
    }
    if (loops < q_.loop_count)
        throw Error(ErrorCode::constraint, "only " + std::to_string(loops) + " of " + std::to_string(q_.loop_count) +
                                               " requested loops could be routed");
}

void Generator::connect_authored_ports() {
    auto hosts_by_distance = [&](Vec3 p) {
        std::vector<std::pair<double, std::size_t>> order;
        for (std::size_t i = 0; i < hosts_.size(); ++i)
            order.push_back({distance(p, out_.caverns[hosts_[i].cavern].center), i});
        std::sort(order.begin(), order.end());
        return order;
    };
    const std::size_t authored_ports = q_.authored_ports.size();
    // Required ports not already used by an authored tunnel.
    for (std::size_t i = 0; i < authored_ports; ++i) {
        if (out_.ports[i].state != PortState::required || port_used_[i]) continue;
        std::string why = "the region has no generated chambers or gateways";
        bool ok = false;
        const auto order = hosts_by_distance(out_.ports[i].position);
        for (std::size_t k = 0; k < order.size() && k < 4 && !ok; ++k) ok = connect_authored(i, order[k].second, why);
        if (!ok) throw Error(ErrorCode::no_route, "required port " + out_.ports[i].id + " cannot be connected: " + why);
    }
    // Authored caverns still outside the network (and not tagged isolated).
    const std::size_t authored_caverns = q_.authored_caverns.size();
    auto isolated_tag = [](const Cavern& c) {
        for (const auto& t : c.tags) {
            std::string k = t;
            for (char& ch : k)
                if (ch >= 'A' && ch <= 'Z') ch = static_cast<char>(ch - 'A' + 'a');
            if (k == "isolated") return true;
        }
        return false;
    };
    for (std::size_t c = 0; c < authored_caverns; ++c) {
        if (isolated_tag(out_.caverns[c])) continue;
        const bool has_network = !hosts_.empty();
        if (has_network && find(c) == find(hosts_[0].cavern)) continue;
        if (!has_network) {
            bool linked = false;
            for (std::size_t p = 0; p < authored_ports && !linked; ++p)
                linked = port_used_[p] && uuid_key(out_.ports[p].cavern_id) == uuid_key(out_.caverns[c].id);
            if (!linked)
                throw Error(ErrorCode::constraint, "authored cavern " + out_.caverns[c].id +
                                                       " is unconnected, not tagged isolated, and the region generates no network");
            continue;
        }
        bool ok = false;
        std::string why = "no unused optional or generated port in its component";
        for (std::size_t p = 0; p < authored_ports && !ok; ++p) {
            const Port& port = out_.ports[p];
            if (port_used_[p] || port.state == PortState::sealed || port.state == PortState::required) continue;
            if (find(cavern_index_.at(uuid_key(port.cavern_id))) != find(c)) continue;
            const auto order = hosts_by_distance(port.position);
            for (std::size_t k = 0; k < order.size() && k < 3 && !ok; ++k) ok = connect_authored(p, order[k].second, why);
        }
        if (!ok)
            throw Error(ErrorCode::no_route, "authored cavern " + out_.caverns[c].id +
                                                 " cannot join the regional network: " + why);
    }
}

void Generator::finish() {
    out_.support_points = lattice_.size();
    try {
        check_region(out_, &budget_);
    } catch (const Error& e) {
        if (e.code == ErrorCode::cancelled || e.code == ErrorCode::budget) throw;
        throw Error(ErrorCode::internal, std::string("generated region failed validation: ") + e.what());
    }
    // Request-level guarantees that Region alone cannot express: reserved boxes
    // and authored clearances are honoured by every generated record.
    const double gap = dv_.gap;
    // Only records created by this run (authored copies may carry a generator marker).
    for (std::size_t ti = q_.authored_tunnels.size(); ti < out_.tunnels.size(); ++ti) {
        const Tunnel& t = out_.tunnels[ti];
        double r = 0;
        for (const auto& p : t.centerline) r = std::max(r, section_radius(p.width, p.height));
        const std::string sa = uuid_key(out_.ports[port_index_.at(uuid_key(t.start_port_id))].cavern_id);
        const std::string sb = uuid_key(out_.ports[port_index_.at(uuid_key(t.end_port_id))].cavern_id);
        for (std::size_t i = 1; i < t.centerline.size(); ++i) {
            const Vec3 a = t.centerline[i - 1].center, b = t.centerline[i].center;
            budget_.charge(1 + q_.reserved.size() + q_.authored_caverns.size());
            for (const Bounds& box : q_.reserved)
                if (segment_box_distance(a, b, box) < r + gap - 1e-6)
                    throw Error(ErrorCode::internal, "generated tunnel " + t.id + " enters a reserved box");
            for (const Cavern& c : q_.authored_caverns) {
                const std::string k = uuid_key(c.id);
                if ((i == 1 && k == sa) || (i + 1 == t.centerline.size() && k == sb)) continue;
                if (segment_box_distance(a, b, c.bounds) < r + c.reserved_clearance - 1e-6)
                    throw Error(ErrorCode::internal, "generated tunnel " + t.id + " enters authored cavern clearance");
            }
        }
    }
    for (std::size_t ci = q_.authored_caverns.size(); ci < out_.caverns.size(); ++ci) {
        const Cavern& c = out_.caverns[ci];
        budget_.charge(1 + q_.reserved.size() + q_.authored_caverns.size());
        for (const Bounds& box : q_.reserved)
            if (box_box_distance(c.bounds, box) < gap - 1e-9)
                throw Error(ErrorCode::internal, "generated cavern " + c.id + " enters a reserved box");
        for (const Cavern& a : q_.authored_caverns)
            if (box_box_distance(c.bounds, a.bounds) < gap + a.reserved_clearance - 1e-9)
                throw Error(ErrorCode::internal, "generated cavern " + c.id + " enters authored cavern clearance");
    }
}

Region Generator::run() {
    validate_request(q_, o_, &budget_);
    budget_.check_cancelled();
    setup();
    copy_authored();
    place_gateways();
    place_chambers();
    build_lattice();
    std::size_t needs = 0;
    for (const Port& p : q_.authored_ports) needs += p.state != PortState::sealed ? 1 : 0;
    max_route_attempts_ = 8ull * hosts_.size() + 4ull * q_.loop_count + 4ull * needs + 16ull;
    connect_network();
    connect_authored_ports();
    finish();
    out_.work_used = budget_.used();
    return std::move(out_);
}

} // namespace
} // namespace detail

Region generate(const Request& request, const Options& options) {
    try {
        detail::Generator generator(request, options);
        return generator.run();
    } catch (const Error&) {
        throw;
    } catch (const std::bad_alloc&) {
        throw Error(ErrorCode::budget, "generation ran out of memory");
    } catch (const std::exception& e) {
        throw Error(ErrorCode::internal, std::string("generation failed: ") + e.what());
    }
}

} // namespace fbs::caves
