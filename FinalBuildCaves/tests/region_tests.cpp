// Regional generator tests (generate / validate / stitch / stable_id).
// Self-contained: no test framework, only the public header fbs/caves.hpp.

#include "fbs/caves.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

using namespace fbs::caves;

namespace {

int failures = 0;
int checks = 0;

#define CHECK(cond)                                                                                  \
    do {                                                                                             \
        ++checks;                                                                                    \
        if (!(cond)) {                                                                               \
            ++failures;                                                                              \
            std::fprintf(stderr, "  FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond);                   \
        }                                                                                            \
    } while (0)

template <class F> ErrorCode error_of(F&& f, std::string* message = nullptr) {
    try {
        f();
    } catch (const Error& e) {
        if (message) *message = e.what();
        return e.code;
    } catch (const std::exception& e) {
        if (message) *message = e.what();
        return static_cast<ErrorCode>(-1);
    }
    return static_cast<ErrorCode>(-2);
}

#define CHECK_ERROR(expr, expected)                                                                   \
    do {                                                                                              \
        std::string msg_;                                                                             \
        const ErrorCode got_ = error_of([&] { expr; }, &msg_);                                        \
        ++checks;                                                                                     \
        if (got_ != (expected)) {                                                                     \
            ++failures;                                                                               \
            std::fprintf(stderr, "  FAIL %s:%d: %s -> code %d (%s), expected %d\n", __FILE__, __LINE__, \
                         #expr, static_cast<int>(got_), msg_.c_str(), static_cast<int>(expected));    \
        }                                                                                             \
    } while (0)

// ----------------------------------------------------------------- dumping
std::string hex(double v) {
    std::uint64_t bits;
    std::memcpy(&bits, &v, 8);
    char b[24];
    std::snprintf(b, sizeof b, "%016" PRIx64, bits);
    return b;
}
std::string dump(Vec3 v) { return hex(v.x) + "," + hex(v.y) + "," + hex(v.z); }
std::string dump(const Polygon& p) {
    std::string s = "[";
    for (auto& v : p) s += hex(v.x) + "," + hex(v.y) + ";";
    return s + "]";
}
std::string dump(const Generation& g) {
    std::string s = std::to_string(g.seed) + "{";
    for (auto& [k, v] : g.parameters) s += k + "=" + v + ";";
    return s + "}";
}
std::string dump(const std::vector<std::string>& t) {
    std::string s = "<";
    for (auto& x : t) s += x + ";";
    return s + ">";
}
std::string dump(const Cavern& c) {
    return c.id + "|" + c.code + "|" + c.name + "|" + dump(c.center) + "|" + dump(c.bounds.minimum) + "|" +
           dump(c.bounds.maximum) + "|" + dump(c.boundary) + "|" + std::to_string(int(c.authoring)) + "|" +
           std::to_string(c.locked) + "|" + hex(c.reserved_clearance) + "|" + dump(c.generation) + "|" + dump(c.tags);
}
std::string dump(const Floor& f) {
    return f.id + "|" + f.cavern_id + "|" + f.code + "|" + f.name + "|" + dump(f.boundary) + "|" + hex(f.base_z) + "|" +
           hex(f.slope_x) + "|" + hex(f.slope_y) + "|" + hex(f.variation_amplitude) + "|" +
           std::to_string(f.variation_seed) + "|" + f.material + "|" + std::to_string(int(f.traversal)) + "|" + dump(f.tags);
}
std::string dump(const Cliff& k) {
    return k.id + "|" + k.cavern_id + "|" + k.from_floor_id + "|" + k.to_floor_id + "|" + k.code + "|" + dump(k.boundary) +
           "|" + std::to_string(int(k.transition)) + "|" + hex(k.height);
}
std::string dump(const Port& p) {
    return p.id + "|" + p.cavern_id + "|" + p.floor_id + "|" + p.code + "|" + p.name + "|" + dump(p.position) + "|" +
           dump(p.facing) + "|" + hex(p.width) + "|" + hex(p.height) + "|" + std::to_string(int(p.state)) + "|" +
           std::to_string(int(p.traversal));
}
std::string dump(const Tunnel& t) {
    std::string s = t.id + "|" + t.start_port_id + "|" + t.end_port_id + "|" + t.code + "|" + t.name + "|";
    for (auto& p : t.centerline) s += dump(p.center) + "," + hex(p.width) + "," + hex(p.height) + ";";
    s += "|" + hex(t.maximum_slope_degrees) + "|" + hex(t.minimum_clearance) + "|" + std::to_string(int(t.traversal)) + "|" +
         std::to_string(t.locked) + "|" + dump(t.generation) + "|" + dump(t.tags) + "|";
    for (auto n : t.source_node_ids) s += std::to_string(n) + ",";
    return s;
}
std::string dump(const Region& r) {
    std::ostringstream o;
    o << r.version << r.algorithm << r.id << r.settings_id << r.world_seed << r.key.x << "," << r.key.y << "," << r.key.z
      << dump(r.size) << dump(r.origin) << "\n";
    for (auto& x : r.caverns) o << dump(x) << "\n";
    for (auto& x : r.floors) o << dump(x) << "\n";
    for (auto& x : r.cliffs) o << dump(x) << "\n";
    for (auto& x : r.ports) o << dump(x) << "\n";
    for (auto& x : r.tunnels) o << dump(x) << "\n";
    for (auto& b : r.boundaries)
        o << b.key << b.port_id << int(b.face) << dump(b.position) << hex(b.width) << hex(b.height) << "\n";
    o << r.support_points << "\n";
    return o.str();
}

// ----------------------------------------------------------------- helpers
Request small(std::int64_t x, std::int64_t y, std::int64_t z, std::uint64_t seed = 42) {
    Request q;
    q.world_seed = seed;
    q.key = {x, y, z};
    q.size = {160, 160, 48};
    q.chamber_count = 5;
    q.loop_count = 1;
    q.chamber_radius_min = 6;
    q.chamber_radius_max = 10;
    q.chamber_height = 8;
    q.tunnel_width = 3;
    q.tunnel_height = 2.5;
    q.sampling_radius = 0.1;
    q.neighbors = 16;
    return q;
}

double dist(Vec3 a, Vec3 b) { return std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z)); }
bool same(Vec3 a, Vec3 b) { return a.x == b.x && a.y == b.y && a.z == b.z; }
double point_box(Vec3 p, const Bounds& b) {
    const double dx = std::max({b.minimum.x - p.x, 0.0, p.x - b.maximum.x});
    const double dy = std::max({b.minimum.y - p.y, 0.0, p.y - b.maximum.y});
    const double dz = std::max({b.minimum.z - p.z, 0.0, p.z - b.maximum.z});
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}
// Independent (sampling based) minimum distance from a tunnel's centreline to a box.
double tunnel_box_distance(const Tunnel& t, const Bounds& b, std::size_t skip_first, std::size_t skip_last) {
    double best = std::numeric_limits<double>::infinity();
    const std::size_t n = t.centerline.size();
    for (std::size_t i = 1 + skip_first; i + skip_last < n; ++i) {
        const Vec3 a = t.centerline[i - 1].center, c = t.centerline[i].center;
        const int steps = std::max(2, static_cast<int>(dist(a, c) / 0.05));
        for (int s = 0; s <= steps; ++s) {
            const double u = double(s) / steps;
            best = std::min(best, point_box({a.x + (c.x - a.x) * u, a.y + (c.y - a.y) * u, a.z + (c.z - a.z) * u}, b));
        }
    }
    return best;
}
double tunnel_radius(const Tunnel& t) {
    double r = 0;
    for (auto& p : t.centerline) r = std::max(r, 0.5 * std::hypot(p.width, p.height));
    return r;
}
const Port* port_by_id(const Region& r, const std::string& id) {
    for (auto& p : r.ports)
        if (p.id == id) return &p;
    return nullptr;
}
std::size_t components(const Region& r) {
    std::map<std::string, std::string> parent;
    for (auto& c : r.caverns) parent[c.id] = c.id;
    std::function<std::string(const std::string&)> find = [&](const std::string& x) {
        return parent[x] == x ? x : parent[x] = find(parent[x]);
    };
    for (auto& t : r.tunnels) parent[find(port_by_id(r, t.start_port_id)->cavern_id)] = find(port_by_id(r, t.end_port_id)->cavern_id);
    std::set<std::string> roots;
    for (auto& c : r.caverns) roots.insert(find(c.id));
    return roots.size();
}
bool cavern_is_generated(const Region& r, const std::string& id) {
    for (auto& c : r.caverns)
        if (c.id == id) return c.generation.parameters.count("fbs.generator") == 1;
    return false;
}
bool uuid_like(const std::string& s) {
    if (s.size() != 36) return false;
    for (std::size_t i = 0; i < 36; ++i) {
        if (i == 8 || i == 13 || i == 18 || i == 23) {
            if (s[i] != '-') return false;
        } else if (!std::isxdigit(static_cast<unsigned char>(s[i]))) {
            return false;
        }
    }
    return true;
}
double slope_deg(Vec3 a, Vec3 b) {
    const double h = std::hypot(b.x - a.x, b.y - a.y);
    return h <= 1e-9 ? 90 : std::atan2(std::fabs(b.z - a.z), h) * 180 / 3.14159265358979323846;
}
const Boundary* boundary_of(const Region& r, Face f) {
    for (auto& b : r.boundaries)
        if (b.face == f) return &b;
    return nullptr;
}

// MOOCoW UnitVector3 normalisation step (unfused products, left-to-right sum).
bool storage_stable(Vec3 v) {
    volatile double xx = v.x * v.x, yy = v.y * v.y, zz = v.z * v.z;
    const double m = std::sqrt(xx + yy + zz);
    return v.x / m == v.x && v.y / m == v.y && v.z / m == v.z;
}

// Checks the structural guarantees of every generated tunnel independently of validate().
void check_generated_geometry(const Region& g) {
    std::set<std::string> used;
    for (const Tunnel& t : g.tunnels) {
        if (!t.generation.parameters.count("fbs.generator")) continue;
        const Port* a = port_by_id(g, t.start_port_id);
        const Port* b = port_by_id(g, t.end_port_id);
        CHECK(a && b);
        if (!a || !b) continue;
        CHECK(same(t.centerline.front().center, a->position));
        CHECK(same(t.centerline.back().center, b->position));
        CHECK(t.centerline.front().height == a->height && t.centerline.back().height == b->height);
        CHECK(used.insert(a->id).second);
        CHECK(used.insert(b->id).second);
        CHECK(!t.source_node_ids.empty());
        for (std::size_t i = 1; i < t.centerline.size(); ++i) {
            const Profile& p = t.centerline[i - 1];
            const Profile& q = t.centerline[i];
            CHECK(slope_deg(p.center, q.center) <= t.maximum_slope_degrees);
            CHECK(slope_deg({p.center.x, p.center.y, p.center.z - p.height / 2}, {q.center.x, q.center.y, q.center.z - q.height / 2}) <=
                  t.maximum_slope_degrees);
            CHECK(q.width > 0 && q.height >= t.minimum_clearance);
        }
        for (const Port* p : {a, b})
            if (p->code.rfind("P", 0) == 0 && cavern_is_generated(g, p->cavern_id)) CHECK(storage_stable(p->facing));
        // Floors sit at the port bottom elevation.
        for (const Port* p : {a, b}) {
            for (auto& f : g.floors)
                if (f.id == p->floor_id) CHECK(std::fabs(p->position.z - p->height / 2 - f.base_z) < 1e-9);
        }
    }
}

// ----------------------------------------------------------------- tests
void test_stable_id() {
    const std::string a = stable_id("alpha"), b = stable_id("alpha"), c = stable_id("alphb"), e = stable_id("");
    CHECK(a == b);
    CHECK(a != c);
    CHECK(uuid_like(a) && uuid_like(e));
    CHECK(a[14] == '8');
    CHECK(std::string("89ab").find(a[19]) != std::string::npos);
    for (auto& ch : a) CHECK(!(ch >= 'A' && ch <= 'F'));
    // Frozen reference value: identities must never change between builds.
    CHECK(stable_id("fbs-caves/reference") == "34e7fde2-d778-8bf1-9a7a-a931e5b24766");
    std::set<std::string> many;
    for (int i = 0; i < 20000; ++i) many.insert(stable_id("n/" + std::to_string(i)));
    CHECK(many.size() == 20000);
}

void test_generate_basic() {
    const Request q = small(0, 0, 0);
    const Region g = generate(q);
    validate(g);
    CHECK(g.caverns.size() == 9);    // 4 gateways + 5 chambers
    CHECK(g.boundaries.size() == 4);
    CHECK(components(g) == 1);
    // Cyclomatic number of the cavern graph equals the requested loops.
    CHECK(g.tunnels.size() + 1 == g.caverns.size() + q.loop_count);
    CHECK(g.support_points > 0 && g.work_used > 0);
    CHECK(uuid_like(g.id) && uuid_like(g.settings_id));
    CHECK(g.origin.x == 0 && g.origin.y == 0 && g.origin.z == 0);
    check_generated_geometry(g);
    std::size_t loops = 0;
    for (auto& t : g.tunnels) {
        loops += std::count(t.tags.begin(), t.tags.end(), "loop");
        CHECK(t.generation.parameters.at("fbs.route.engine") == "KarstNSim::run_simulation_memory");
    }
    CHECK(loops == q.loop_count);
    for (auto& c : g.cliffs) {
        CHECK(c.transition == Transition::climbable_cliff || c.transition == Transition::stairs);
        CHECK(c.height > 0 && !c.to_floor_id.empty());
    }
    for (auto& c : g.caverns) {
        CHECK(c.boundary.size() >= 10);
        CHECK(c.generation.parameters.count("fbs.world_seed") == 1);
        CHECK(c.code.rfind("R0_0_0-", 0) == 0);
    }
    std::printf("  basic: %zu caverns, %zu floors, %zu cliffs, %zu ports, %zu tunnels, %" PRIu64 " support points\n",
                g.caverns.size(), g.floors.size(), g.cliffs.size(), g.ports.size(), g.tunnels.size(), g.support_points);
}

void test_determinism_and_order() {
    const std::string a1 = dump(generate(small(0, 0, 0)));
    const std::string b1 = dump(generate(small(1, 0, 0)));
    const std::string b2 = dump(generate(small(1, 0, 0)));
    const std::string a2 = dump(generate(small(0, 0, 0)));
    CHECK(a1 == a2);
    CHECK(b1 == b2);
    CHECK(a1 != b1);
    const std::string other_seed = dump(generate(small(0, 0, 0, 43)));
    CHECK(other_seed != a1);
    Request styled = small(0, 0, 0);
    styled.style_id = 2;
    CHECK(generate(styled).settings_id != generate(small(0, 0, 0)).settings_id);
}

void test_neighbors_and_stitch() {
    const Region a = generate(small(0, 0, 0));
    const Region b = generate(small(1, 0, 0));
    const Region c = generate(small(0, 1, 0));
    const Boundary* ab = boundary_of(a, Face::positive_x);
    const Boundary* ba = boundary_of(b, Face::negative_x);
    CHECK(ab && ba);
    if (!ab || !ba) return;
    CHECK(ab->key == ba->key && ab->width == ba->width && ab->height == ba->height);
    CHECK(ab->position.x == a.size.x && ba->position.x == 0);
    CHECK(ab->position.y == ba->position.y && ab->position.z == ba->position.z);
    CHECK(ab->port_id != ba->port_id);

    const auto t1 = stitch(a, b);
    const auto t2 = stitch(b, a);
    CHECK(t1.size() == 1 && t2.size() == 1);
    if (t1.size() == 1 && t2.size() == 1) {
        const Tunnel& t = t1[0];
        CHECK(t.id == t2[0].id);
        CHECK(t.start_port_id == ab->port_id && t.end_port_id == ba->port_id);
        CHECK(t2[0].start_port_id == ba->port_id && t2[0].end_port_id == ab->port_id);
        CHECK(same(t.centerline.front().center, port_by_id(a, ab->port_id)->position));
        CHECK(same(t.centerline[1].center, ab->position));
        const Vec3 bp = port_by_id(b, ba->port_id)->position;
        CHECK(same(t.centerline.back().center, Vec3{bp.x + a.size.x, bp.y, bp.z}));
        CHECK(same(t2[0].centerline.front().center, bp));
        CHECK(same(t2[0].centerline[1].center, ba->position));
        for (auto& p : t.centerline) CHECK(p.width == ab->width && p.height == ab->height);
        CHECK(stitch(a, b)[0].centerline.size() == t.centerline.size());
        CHECK(dump(stitch(a, b)[0]) == dump(t));
    }
    const auto ty = stitch(a, c);
    CHECK(ty.size() == 1 && same(ty[0].centerline[1].center, boundary_of(a, Face::positive_y)->position));

    // Rejections.
    const Region far = generate(small(2, 0, 0));
    CHECK_ERROR(stitch(a, far), ErrorCode::incompatible_boundary);
    const Region diagonal = generate(small(1, 1, 0));
    CHECK_ERROR(stitch(a, diagonal), ErrorCode::incompatible_boundary);
    CHECK_ERROR(stitch(a, a), ErrorCode::incompatible_boundary);
    const Region other_seed = generate(small(1, 0, 0, 7));
    CHECK_ERROR(stitch(a, other_seed), ErrorCode::incompatible_boundary);
    Request wider = small(1, 0, 0);
    wider.tunnel_width = 3.25;
    const Region other_settings = generate(wider);
    CHECK_ERROR(stitch(a, other_settings), ErrorCode::incompatible_boundary);
    Region tampered = b;
    for (auto& bd : tampered.boundaries)
        if (bd.face == Face::negative_x) bd.position.y += 0.5;
    CHECK_ERROR(stitch(a, tampered), ErrorCode::incompatible_boundary);
    Region no_face = b;
    no_face.boundaries.clear();
    CHECK_ERROR(stitch(a, no_face), ErrorCode::incompatible_boundary);
}

void test_large_and_negative_keys() {
    const Region n1 = generate(small(-7, -3, 0));
    const Region n2 = generate(small(-6, -3, 0));
    validate(n1);
    CHECK(n1.origin.x == -7 * 160.0 && n1.origin.y == -3 * 160.0);
    CHECK(stitch(n1, n2).size() == 1);
    CHECK(stitch(n2, n1).size() == 1);

    const std::int64_t max = std::numeric_limits<std::int64_t>::max();
    const std::int64_t min = std::numeric_limits<std::int64_t>::min();
    const Region big = generate(small(max - 1, min + 1, 0));
    validate(big);
    CHECK(std::isfinite(big.origin.x) && std::isfinite(big.origin.y));
    // No neighbour exists beyond the signed range: requesting that face is invalid.
    CHECK_ERROR(generate(small(max, 0, 0)), ErrorCode::invalid_input);
    CHECK_ERROR(generate(small(0, min, 0)), ErrorCode::invalid_input);
    Request edge = small(max, 0, 0);
    edge.face_mask = 1 | 4 | 8; // nx, ny, py
    const Region e = generate(edge);
    const Region inner = generate(small(max - 1, 0, 0));
    CHECK(stitch(inner, e).size() == 1);
    CHECK(stitch(e, inner).size() == 1);
}

void test_concurrency() {
    const int n = 8;
    std::vector<std::string> sequential(n), parallel(n);
    for (int i = 0; i < n; ++i) sequential[i] = dump(generate(small(i % 4, i / 4, 0)));
    std::vector<std::thread> threads;
    for (int i = n - 1; i >= 0; --i)
        threads.emplace_back([&, i] { parallel[i] = dump(generate(small(i % 4, i / 4, 0))); });
    for (auto& t : threads) t.join();
    for (int i = 0; i < n; ++i) CHECK(sequential[i] == parallel[i]);
    // Same region on several threads at once.
    std::vector<std::string> same_region(4);
    threads.clear();
    for (int i = 0; i < 4; ++i) threads.emplace_back([&, i] { same_region[i] = dump(generate(small(3, 3, 0))); });
    for (auto& t : threads) t.join();
    for (int i = 1; i < 4; ++i) CHECK(same_region[i] == same_region[0]);
}

// Authored fixture: a locked exact vault with a required and a sealed port, a
// second cavern joined to it by an authored tunnel, an isolated authored cavern
// and a reserved box.
Request authored_request() {
    Request q = small(0, 0, 0);
    Cavern v;
    v.id = "0f1e2d3c-4b5a-4978-8a6b-5c4d3e2f1a00";
    v.code = "VAULT";
    v.name = "Vault";
    v.center = {70, 70, 15};
    v.bounds = {{60, 60, 10}, {80, 80, 20}};
    v.boundary = {{60, 60}, {80, 60}, {80, 80}, {60, 80}};
    v.authoring = Authoring::exact;
    v.locked = true;
    v.reserved_clearance = 2;
    v.generation.seed = 7;
    v.generation.parameters = {{"author", "max"}, {"note", "hand placed"}};
    v.tags = {"vault", "Story"};
    Cavern w = v;
    w.id = "0F1E2D3C-4B5A-4978-8A6B-5C4D3E2F1A01"; // upper case must be preserved verbatim
    w.code = "ANNEX";
    w.name = "Annex";
    w.center = {70, 100, 15};
    w.bounds = {{60, 92, 10}, {80, 112, 20}};
    w.boundary = {{60, 92}, {80, 92}, {80, 112}, {60, 112}};
    w.locked = false;
    w.authoring = Authoring::constrained_procedural;
    w.tags = {};
    Cavern iso = v;
    iso.id = "0f1e2d3c-4b5a-4978-8a6b-5c4d3e2f1a02";
    iso.code = "SEALED-ROOM";
    iso.name = "Sealed room";
    iso.center = {130, 130, 30};
    iso.bounds = {{125, 125, 26}, {135, 135, 34}};
    iso.boundary = {};
    iso.authoring = Authoring::constrained_procedural;
    iso.tags = {"Isolated"};
    q.authored_caverns = {v, w, iso};

    Floor f;
    f.id = "1f1e2d3c-4b5a-4978-8a6b-5c4d3e2f1a00";
    f.cavern_id = v.id;
    f.code = "F";
    f.name = "Vault floor";
    f.boundary = v.boundary;
    f.base_z = 10;
    f.slope_x = 0.01;
    f.variation_amplitude = 0.2;
    f.variation_seed = 99;
    f.material = "granite";
    f.traversal = Traversal::difficult;
    f.tags = {"authored"};
    Floor g = f;
    g.id = "1f1e2d3c-4b5a-4978-8a6b-5c4d3e2f1a01";
    g.cavern_id = w.id;
    g.boundary = w.boundary;
    g.name = "Annex floor";
    q.authored_floors = {f, g};

    Cliff k;
    k.id = "2f1e2d3c-4b5a-4978-8a6b-5c4d3e2f1a00";
    k.cavern_id = v.id;
    k.from_floor_id = f.id;
    k.to_floor_id = "";
    k.code = "EDGE";
    k.boundary = {{62, 62}, {62, 78}};
    k.transition = Transition::impassable_cliff;
    k.height = 4;
    q.authored_cliffs = {k};

    Port required;
    required.id = "3f1e2d3c-4b5a-4978-8a6b-5c4d3e2f1a00";
    required.cavern_id = v.id;
    required.floor_id = f.id;
    required.code = "EAST";
    required.name = "East door";
    required.position = {79.5, 70, 11.25};
    required.facing = {1, 0, 0};
    required.width = 3;
    required.height = 2.5;
    required.state = PortState::required;
    required.traversal = Traversal::walkable;
    Port sealed = required;
    sealed.id = "3f1e2d3c-4b5a-4978-8a6b-5c4d3e2f1a01";
    sealed.code = "WEST";
    sealed.name = "West door";
    sealed.position = {60.5, 70, 11.25};
    sealed.facing = {-1, 0, 0};
    sealed.state = PortState::sealed;
    Port north = required;
    north.id = "3f1e2d3c-4b5a-4978-8a6b-5c4d3e2f1a02";
    north.code = "NORTH";
    north.name = "North door";
    north.position = {70, 79.5, 11.25};
    north.facing = {0, 1, 0};
    north.state = PortState::optional;
    Port annex = north;
    annex.id = "3f1e2d3c-4b5a-4978-8a6b-5c4d3e2f1a03";
    annex.cavern_id = w.id;
    annex.floor_id = g.id;
    annex.code = "SOUTH";
    annex.name = "Annex door";
    annex.position = {70, 92.5, 11.25};
    annex.facing = {0, -1, 0};
    q.authored_ports = {required, sealed, north, annex};

    Tunnel t;
    t.id = "4f1e2d3c-4b5a-4978-8a6b-5c4d3e2f1a00";
    t.start_port_id = north.id;
    t.end_port_id = annex.id;
    t.code = "LINK";
    t.name = "Vault link";
    t.centerline = {{{70, 79.5, 11.25}, 3, 2.5}, {{70, 86, 11.25}, 3.5, 2.75}, {{70, 92.5, 11.25}, 3, 2.5}};
    t.maximum_slope_degrees = 45;
    t.minimum_clearance = 2.5;
    t.traversal = Traversal::restricted;
    t.locked = true;
    t.generation.seed = 11;
    t.generation.parameters = {{"author", "max"}};
    t.tags = {"authored"};
    t.source_node_ids = {5, 6};
    q.authored_tunnels = {t};
    q.reserved = {{{100, 32, 0}, {124, 54, 48}}};
    return q;
}

void test_authored_preservation() {
    const Request q = authored_request();
    const Region g = generate(q);
    validate(g);
    check_generated_geometry(g);
    // Exact copies, in input order, before generated records.
    for (std::size_t i = 0; i < q.authored_caverns.size(); ++i) CHECK(dump(g.caverns[i]) == dump(q.authored_caverns[i]));
    for (std::size_t i = 0; i < q.authored_floors.size(); ++i) CHECK(dump(g.floors[i]) == dump(q.authored_floors[i]));
    for (std::size_t i = 0; i < q.authored_cliffs.size(); ++i) CHECK(dump(g.cliffs[i]) == dump(q.authored_cliffs[i]));
    for (std::size_t i = 0; i < q.authored_ports.size(); ++i) CHECK(dump(g.ports[i]) == dump(q.authored_ports[i]));
    for (std::size_t i = 0; i < q.authored_tunnels.size(); ++i) CHECK(dump(g.tunnels[i]) == dump(q.authored_tunnels[i]));
    // Required port connected exactly once by a generated tunnel; sealed port unused.
    int required_uses = 0, sealed_uses = 0;
    for (auto& t : g.tunnels) {
        for (auto* id : {&t.start_port_id, &t.end_port_id}) {
            required_uses += *id == q.authored_ports[0].id;
            sealed_uses += *id == q.authored_ports[1].id;
        }
    }
    CHECK(required_uses == 1);
    CHECK(sealed_uses == 0);
    // Generated geometry keeps clear of reserved and authored volumes (independent sampling check).
    for (auto& t : g.tunnels) {
        if (!t.generation.parameters.count("fbs.generator")) continue;
        const double r = tunnel_radius(t);
        CHECK(tunnel_box_distance(t, q.reserved[0], 0, 0) >= r);
        for (auto& c : q.authored_caverns) {
            const bool start = port_by_id(g, t.start_port_id)->cavern_id == c.id;
            const bool end = port_by_id(g, t.end_port_id)->cavern_id == c.id;
            CHECK(tunnel_box_distance(t, c.bounds, start ? 1 : 0, end ? 1 : 0) >= r + c.reserved_clearance);
        }
        // The authored passage is protected like any other volume.
        const Bounds link{{70, 79.5, 11.25}, {70, 92.5, 11.25}};
        CHECK(tunnel_box_distance(t, link, 0, 0) >= r + 0.5 * std::hypot(3.5, 2.75));
    }
    for (auto& c : g.caverns) {
        if (!c.generation.parameters.count("fbs.generator")) continue;
        CHECK(point_box(c.bounds.minimum, q.reserved[0]) > 0 || point_box(c.bounds.maximum, q.reserved[0]) > 0);
        const Bounds& rb = q.reserved[0];
        const bool overlap = c.bounds.minimum.x <= rb.maximum.x && c.bounds.maximum.x >= rb.minimum.x &&
                             c.bounds.minimum.y <= rb.maximum.y && c.bounds.maximum.y >= rb.minimum.y;
        CHECK(!overlap);
    }
    CHECK(components(g) == 2); // everything except the isolated room is one network
    // The authored tunnel is protected: its segment keeps generated passages at a distance.
    // Deterministic with authored content too.
    CHECK(dump(generate(q)) == dump(g));
}

void test_authored_failures() {
    {
        Request q = authored_request();
        q.authored_ports[0].facing = {0, 0, 1}; // required port facing straight up
        CHECK_ERROR(generate(q), ErrorCode::no_route);
    }
    {
        Request q = authored_request();
        q.authored_tunnels[0].end_port_id = q.authored_ports[1].id; // uses sealed port
        q.authored_tunnels[0].start_port_id = q.authored_ports[2].id;
        CHECK_ERROR(generate(q), ErrorCode::invalid_input);
    }
    {
        Request q = authored_request();
        q.authored_caverns[0].id = "not-a-uuid";
        CHECK_ERROR(generate(q), ErrorCode::invalid_input);
    }
    {
        Request q = authored_request();
        q.authored_caverns[0].id = "00000000-0000-0000-0000-000000000000";
        CHECK_ERROR(generate(q), ErrorCode::invalid_input);
    }
    {
        Request q = authored_request();
        q.authored_floors[0].boundary = {{60, 60}, {80, 80}, {80, 60}, {60, 80}}; // bow tie
        CHECK_ERROR(generate(q), ErrorCode::invalid_input);
    }
    {
        Request q = authored_request();
        q.authored_ports[0].cavern_id = "9f1e2d3c-4b5a-4978-8a6b-5c4d3e2f1a00"; // dangling reference
        CHECK_ERROR(generate(q), ErrorCode::invalid_input);
    }
    {
        Request q = authored_request();
        q.authored_ports[3].id = q.authored_ports[2].id; // duplicate id
        CHECK_ERROR(generate(q), ErrorCode::invalid_input);
    }
    {
        Request q = authored_request();
        q.authored_caverns[1].bounds.maximum.x = 170; // outside region
        CHECK_ERROR(generate(q), ErrorCode::invalid_input);
    }
    {
        Request q = authored_request();
        q.authored_tunnels[0].centerline[1].center.z = 40; // steeper than its own limit
        CHECK_ERROR(generate(q), ErrorCode::invalid_input);
    }
    {
        Request q = authored_request();
        q.reserved.push_back({{0, 0, 0}, {160, 160, 48}}); // everything reserved
        CHECK_ERROR(generate(q), ErrorCode::constraint);
    }
    {
        // Required port whose access stub would have to leave the region.
        Request q = authored_request();
        Cavern edge = q.authored_caverns[2];
        edge.id = "0f1e2d3c-4b5a-4978-8a6b-5c4d3e2f1a09";
        edge.code = "EDGE";
        edge.center = {154, 75, 15};
        edge.bounds = {{150, 70, 10}, {158, 80, 20}};
        edge.tags = {};
        Port p = q.authored_ports[0];
        p.id = "3f1e2d3c-4b5a-4978-8a6b-5c4d3e2f1a09";
        p.cavern_id = edge.id;
        p.floor_id = "";
        p.position = {157.5, 75, 11.25};
        q.authored_caverns.push_back(edge);
        q.authored_ports.push_back(p);
        q.face_mask = 1;
        CHECK_ERROR(generate(q), ErrorCode::no_route);
    }
    {
        // A required port on a cavern with an absurd reservation cannot be reached.
        Request q = authored_request();
        q.face_mask = 0;
        q.chamber_count = 1;
        q.loop_count = 0;
        q.authored_caverns[0].reserved_clearance = 70;
        CHECK_ERROR(generate(q), ErrorCode::no_route);
    }
}

void test_authored_text_rules() {
    {
        Request q = authored_request();
        q.authored_caverns[0].code = " VAULT";
        CHECK_ERROR(generate(q), ErrorCode::invalid_input);
    }
    {
        Request q = authored_request();
        q.authored_caverns[0].name = "Vault\u00a0"; // trailing no-break space would be trimmed by MOOCoW
        CHECK_ERROR(generate(q), ErrorCode::invalid_input);
    }
    {
        Request q = authored_request();
        q.authored_ports[0].name = "East\ndoor";
        CHECK_ERROR(generate(q), ErrorCode::invalid_input);
    }
    {
        Request q = authored_request();
        q.authored_caverns[0].tags = {"vault", "VAULT"};
        CHECK_ERROR(generate(q), ErrorCode::invalid_input);
    }
    {
        Request q = authored_request();
        q.authored_caverns[0].generation.parameters["parameters"] = "{}";
        CHECK_ERROR(generate(q), ErrorCode::invalid_input);
    }
    {
        // Floor and port share one cavern code scope in MOOCoW.
        Request q = authored_request();
        q.authored_floors[0].code = "EAST";
        CHECK_ERROR(generate(q), ErrorCode::invalid_input);
    }
    {
        // Cavern and tunnel share the map code scope.
        Request q = authored_request();
        q.authored_tunnels[0].code = "vault";
        CHECK_ERROR(generate(q), ErrorCode::invalid_input);
    }
    {
        // An authored code equal to a generated code: the generated record yields.
        Request q = authored_request();
        q.authored_caverns[1].code = "r0_0_0-c01";
        const Region g = generate(q);
        validate(g);
        CHECK(g.caverns[1].code == "r0_0_0-c01");
        bool renamed = false;
        for (auto& c : g.caverns) renamed = renamed || c.code == "R0_0_0-C01~2";
        CHECK(renamed);
    }
}

// Previously generated records fed back as locked authored content keep their
// identity and never collide with the regenerated remainder.
void test_locked_generated_reuse() {
    const Request base = small(0, 0, 0);
    const Region first = generate(base);
    const Cavern* chamber = nullptr;
    for (auto& c : first.caverns)
        if (c.code == "R0_0_0-C02") chamber = &c;
    CHECK(chamber != nullptr);
    if (!chamber) return;
    Request q = base;
    Cavern locked = *chamber;
    locked.locked = true;
    q.authored_caverns.push_back(locked);
    for (auto& f : first.floors)
        if (f.cavern_id == locked.id) q.authored_floors.push_back(f);
    for (auto& k : first.cliffs)
        if (k.cavern_id == locked.id) q.authored_cliffs.push_back(k);
    for (auto& p : first.ports)
        if (p.cavern_id == locked.id) q.authored_ports.push_back(p);
    const Region g = generate(q);
    validate(g);
    CHECK(dump(g.caverns[0]) == dump(locked));
    for (std::size_t i = 0; i < q.authored_ports.size(); ++i) CHECK(dump(g.ports[i]) == dump(q.authored_ports[i]));
    std::set<std::string> ids;
    for (auto& c : g.caverns) CHECK(ids.insert(c.id).second);
    for (auto& t : g.tunnels) CHECK(ids.insert(t.id).second);
    for (auto& p : g.ports) CHECK(ids.insert(p.id).second);
    CHECK(components(g) == 1);
    check_generated_geometry(g);
}

void test_loops() {
    Request q = small(0, 0, 0);
    q.face_mask = 0;
    q.chamber_count = 7;
    q.loop_count = 3;
    const Region g = generate(q);
    validate(g);
    CHECK(components(g) == 1);
    CHECK(g.tunnels.size() + 1 == g.caverns.size() + 3);
    q.loop_count = 0;
    const Region tree = generate(q);
    CHECK(tree.tunnels.size() + 1 == tree.caverns.size());
}

void test_impossible_constraints() {
    Request q = small(0, 0, 0);
    q.chamber_count = 200;
    CHECK_ERROR(generate(q), ErrorCode::constraint);
    q = small(0, 0, 0);
    q.face_mask = 0;
    q.chamber_count = 2;
    q.loop_count = 5;
    CHECK_ERROR(generate(q), ErrorCode::constraint);
    q = small(0, 0, 0);
    q.chamber_height = 47;
    CHECK_ERROR(generate(q), ErrorCode::constraint);
    q = small(0, 0, 0);
    q.tunnel_width = std::nan("");
    CHECK_ERROR(generate(q), ErrorCode::invalid_input);
    q = small(0, 0, 0);
    q.face_mask = 64;
    CHECK_ERROR(generate(q), ErrorCode::invalid_input);
    q = small(0, 0, 0);
    q.version = schema_version + 1;
    CHECK_ERROR(generate(q), ErrorCode::incompatible_version);
    q = small(0, 0, 0);
    q.size = {0, 160, 48};
    CHECK_ERROR(generate(q), ErrorCode::invalid_input);
    q = small(0, 0, 0);
    q.size = {1e300, 160, 48};
    CHECK_ERROR(generate(q), ErrorCode::invalid_input);
    q = small(0, 0, 0);
    q.maximum_slope_degrees = 90;
    CHECK_ERROR(generate(q), ErrorCode::invalid_input);
    q = small(0, 0, 0);
    q.chamber_radius_min = 20;
    q.chamber_radius_max = 10;
    CHECK_ERROR(generate(q), ErrorCode::invalid_input);
    q = small(0, 0, 0);
    q.neighbors = 1000;
    CHECK_ERROR(generate(q), ErrorCode::invalid_input);
    // Vertical face with a slope limit too shallow for a ramp that fits.
    q = small(0, 0, 0);
    q.face_mask = 32;
    q.maximum_slope_degrees = 2;
    CHECK_ERROR(generate(q), ErrorCode::constraint);
    // Lattice larger than the point budget.
    q = small(0, 0, 0);
    q.sampling_radius = 0.001;
    CHECK_ERROR(generate(q), ErrorCode::budget);
}

void test_cancellation_and_budget() {
    Options cancel_now;
    cancel_now.cancelled = [] { return true; };
    CHECK_ERROR(generate(small(0, 0, 0), cancel_now), ErrorCode::cancelled);
    std::atomic<int> calls{0};
    Options cancel_later;
    cancel_later.cancelled = [&] { return ++calls > 40; };
    CHECK_ERROR(generate(small(0, 0, 0), cancel_later), ErrorCode::cancelled);
    CHECK(calls > 40);
    Options tiny;
    tiny.work_limit = 5000;
    CHECK_ERROR(generate(small(0, 0, 0), tiny), ErrorCode::budget);
    Options few_points;
    few_points.max_points = 100;
    CHECK_ERROR(generate(small(0, 0, 0), few_points), ErrorCode::budget);
    Options zero;
    zero.work_limit = 0;
    CHECK_ERROR(generate(small(0, 0, 0), zero), ErrorCode::invalid_input);
    // A generous budget succeeds and reports its usage.
    const Region g = generate(small(0, 0, 0));
    Options exact;
    exact.work_limit = g.work_used;
    CHECK(generate(small(0, 0, 0), exact).work_used == g.work_used);
    // One unit less must fail: native route work is charged to the same budget.
    Options short_by_one;
    short_by_one.work_limit = g.work_used - 1;
    CHECK_ERROR(generate(small(0, 0, 0), short_by_one), ErrorCode::budget);
}

void test_vertical_faces() {
    Request lower = small(0, 0, 0);
    lower.face_mask = 15 | 32;
    Request upper = small(0, 0, 1);
    upper.face_mask = 15 | 16;
    const Region a = generate(lower);
    const Region b = generate(upper);
    validate(a);
    validate(b);
    const auto t = stitch(a, b);
    CHECK(t.size() == 1);
    if (t.size() != 1) return;
    CHECK(t[0].centerline[1].center.z == a.size.z);
    CHECK(t[0].centerline.back().center.z > a.size.z);
    CHECK(t[0].centerline.front().center.z < a.size.z);
    for (std::size_t i = 1; i < t[0].centerline.size(); ++i)
        CHECK(slope_deg(t[0].centerline[i - 1].center, t[0].centerline[i].center) <= lower.maximum_slope_degrees);
    const auto r = stitch(b, a);
    CHECK(r.size() == 1 && r[0].centerline[1].center.z == 0 && r[0].id == t[0].id);
    // X neighbours still stitch when Z faces are enabled.
    Request side = small(1, 0, 0);
    side.face_mask = 15 | 32;
    CHECK(stitch(a, generate(side)).size() == 1);
    // Mismatched Z enablement: the neighbour lacks the face.
    CHECK_ERROR(stitch(a, generate(small(0, 0, 1))), ErrorCode::incompatible_boundary);
}

void test_validate_detects_tampering() {
    const Region g = generate(small(0, 0, 0));
    std::size_t gen = 0;
    while (gen < g.tunnels.size() && !g.tunnels[gen].generation.parameters.count("fbs.generator")) ++gen;
    {
        Region r = g;
        r.tunnels[gen].centerline.front().center.x += 1e-3;
        CHECK_ERROR(validate(r), ErrorCode::invalid_input);
    }
    {
        Region r = g;
        r.ports[1].id = r.ports[0].id;
        CHECK_ERROR(validate(r), ErrorCode::invalid_input);
    }
    {
        Region r = g;
        r.caverns[1].code = r.caverns[0].code;
        CHECK_ERROR(validate(r), ErrorCode::invalid_input);
    }
    {
        Region r = g;
        r.floors[0].cavern_id = "9f1e2d3c-4b5a-4978-8a6b-5c4d3e2f1a00";
        CHECK_ERROR(validate(r), ErrorCode::invalid_input);
    }
    {
        Region r = g;
        for (auto& p : r.ports)
            if (p.id == r.tunnels[gen].start_port_id) p.state = PortState::sealed;
        CHECK_ERROR(validate(r), ErrorCode::invalid_input);
    }
    {
        // Disconnect a chamber: drop every tunnel touching it.
        Region r = g;
        const std::string victim = r.caverns[5].id;
        r.tunnels.erase(std::remove_if(r.tunnels.begin(), r.tunnels.end(),
                                       [&](const Tunnel& t) {
                                           return port_by_id(g, t.start_port_id)->cavern_id == victim ||
                                                  port_by_id(g, t.end_port_id)->cavern_id == victim;
                                       }),
                        r.tunnels.end());
        CHECK_ERROR(validate(r), ErrorCode::invalid_input);
    }
    {
        // Two generated passages forced to overlap.
        Region r = g;
        Tunnel copy = r.tunnels[gen];
        copy.id = "5f1e2d3c-4b5a-4978-8a6b-5c4d3e2f1a00";
        copy.code = "COPY";
        Port pa = *port_by_id(g, copy.start_port_id), pb = *port_by_id(g, copy.end_port_id);
        pa.id = "6f1e2d3c-4b5a-4978-8a6b-5c4d3e2f1a00";
        pa.code = "PX";
        pb.id = "6f1e2d3c-4b5a-4978-8a6b-5c4d3e2f1a01";
        pb.code = "PY";
        copy.start_port_id = pa.id;
        copy.end_port_id = pb.id;
        r.ports.push_back(pa);
        r.ports.push_back(pb);
        r.tunnels.push_back(copy);
        CHECK_ERROR(validate(r), ErrorCode::invalid_input);
    }
    {
        Region r = g;
        r.boundaries[0].position.z += 0.25;
        CHECK_ERROR(validate(r), ErrorCode::invalid_input);
    }
    {
        Region r = g;
        r.tunnels[gen].centerline[1].center.z += 30; // steep
        CHECK_ERROR(validate(r), ErrorCode::invalid_input);
    }
    {
        Region r = g;
        r.settings_id = stable_id("other");
        CHECK_ERROR(validate(r), ErrorCode::invalid_input);
    }
    {
        Region r = g;
        r.version = 99;
        CHECK_ERROR(validate(r), ErrorCode::incompatible_version);
    }
}

// Passages that cross in plan view keep distinct identities and 3D separation.
void test_crossings() {
    std::size_t crossings = 0;
    for (std::int64_t k = 0; k < 4; ++k) {
        Request q = small(k, 5, 0);
        q.chamber_count = 8;
        q.loop_count = 3;
        const Region g = generate(q);
        validate(g);
        for (std::size_t i = 0; i < g.tunnels.size(); ++i)
            for (std::size_t j = i + 1; j < g.tunnels.size(); ++j) {
                const Tunnel& a = g.tunnels[i];
                const Tunnel& b = g.tunnels[j];
                for (std::size_t s = 1; s < a.centerline.size(); ++s)
                    for (std::size_t t = 1; t < b.centerline.size(); ++t) {
                        const Vec3 p = a.centerline[s - 1].center, p2 = a.centerline[s].center;
                        const Vec3 r = b.centerline[t - 1].center, r2 = b.centerline[t].center;
                        auto orient = [](Vec3 a1, Vec3 b1, Vec3 c1) {
                            return (b1.x - a1.x) * (c1.y - a1.y) - (b1.y - a1.y) * (c1.x - a1.x);
                        };
                        if (orient(p, p2, r) * orient(p, p2, r2) < 0 && orient(r, r2, p) * orient(r, r2, p2) < 0) {
                            ++crossings;
                            CHECK(a.start_port_id != b.start_port_id && a.end_port_id != b.end_port_id);
                        }
                    }
            }
    }
    std::printf("  plan-view crossings observed (all validated as 3D-separated): %zu\n", crossings);
}

// Review regressions: strict UTF-8, bounded validation of many overlapping
// authored passages, and route windows when the lattice margin exceeds 0.7 s.
void test_review_regressions() {
    {
        Request q = authored_request();
        q.authored_caverns[0].name = "Vault \xff";
        CHECK_ERROR(generate(q), ErrorCode::invalid_input);
    }
    {
        Request q = authored_request();
        q.authored_caverns[0].generation.parameters["note"] = "\xc0\x80"; // overlong, even in an optional value
        CHECK_ERROR(generate(q), ErrorCode::invalid_input);
    }
    {
        Region r = generate(small(0, 0, 0));
        r.caverns[0].tags.push_back("\xed\xa0\x80"); // UTF-8 encoded surrogate
        CHECK_ERROR(validate(r), ErrorCode::invalid_input);
    }
    {
        // 20000 authored passages overlapping each other inside one isolated
        // authored cavern: authored/authored pairs are exempt and must not be
        // enumerated, so public validate() stays fast and inside its budget.
        Region r = generate(small(0, 0, 0));
        Cavern c;
        c.id = stable_id("test/overlap/cavern");
        c.code = "OVERLAP";
        c.name = "Overlap";
        c.center = {5, 5, 2};
        c.bounds = {{2, 2, 0.5}, {8, 8, 4}};
        c.authoring = Authoring::constrained_procedural;
        c.tags = {"isolated"};
        Port a;
        a.id = stable_id("test/overlap/a");
        a.cavern_id = c.id;
        a.code = "A";
        a.name = "A";
        a.position = {3, 5, 2};
        a.facing = {-1, 0, 0};
        a.state = PortState::optional;
        Port b = a;
        b.id = stable_id("test/overlap/b");
        b.code = "B";
        b.name = "B";
        b.position = {7, 5, 2};
        b.facing = {1, 0, 0};
        r.caverns.push_back(c);
        r.ports.push_back(a);
        r.ports.push_back(b);
        bool clear_of_generated = true;
        for (auto& gc : r.caverns)
            if (gc.id != c.id && gc.bounds.minimum.x <= 12 && gc.bounds.minimum.y <= 12) clear_of_generated = false;
        CHECK(clear_of_generated);
        for (int i = 0; i < 20000; ++i) {
            Tunnel t;
            t.id = stable_id("test/overlap/t/" + std::to_string(i));
            t.start_port_id = a.id;
            t.end_port_id = b.id;
            t.code = "OV" + std::to_string(i);
            t.name = "Overlap " + std::to_string(i);
            t.centerline = {{a.position, 2, 2}, {b.position, 2, 2}};
            t.maximum_slope_degrees = 45;
            t.minimum_clearance = 2;
            r.tunnels.push_back(t);
        }
        const auto t0 = std::chrono::steady_clock::now();
        validate(r);
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::printf("  validate with 20000 overlapping authored passages: %.2f s\n", seconds);
        CHECK(seconds < 10.0);
    }
    {
        // Wide passages: lattice margin (r + gap/2 ~ 5 m) exceeds 0.7 * spacing (4.8 m).
        Request q = small(0, 0, 0);
        q.tunnel_width = 6;
        q.tunnel_height = 3;
        q.chamber_radius_min = 8;
        q.chamber_radius_max = 12;
        q.chamber_count = 3;
        const Region g = generate(q);
        validate(g);
        check_generated_geometry(g);
        CHECK(components(g) == 1);
    }
}

void run(const char* name, void (*fn)()) {
    const int before = failures;
    const auto t0 = std::chrono::steady_clock::now();
    std::printf("[ RUN  ] %s\n", name);
    std::fflush(stdout);
    try {
        fn();
    } catch (const Error& e) {
        ++failures;
        std::fprintf(stderr, "  FAIL unexpected fbs::caves::Error code %d: %s\n", int(e.code), e.what());
    } catch (const std::exception& e) {
        ++failures;
        std::fprintf(stderr, "  FAIL unexpected exception: %s\n", e.what());
    }
    const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    std::printf("[ %s ] %s (%.0f ms)\n", failures == before ? " OK " : "FAIL", name, ms);
}

} // namespace

int main() {
    run("stable_id", test_stable_id);
    run("generate_basic", test_generate_basic);
    run("determinism_and_order", test_determinism_and_order);
    run("neighbors_and_stitch", test_neighbors_and_stitch);
    run("large_and_negative_keys", test_large_and_negative_keys);
    run("concurrency", test_concurrency);
    run("authored_preservation", test_authored_preservation);
    run("authored_failures", test_authored_failures);
    run("authored_text_rules", test_authored_text_rules);
    run("locked_generated_reuse", test_locked_generated_reuse);
    run("loops", test_loops);
    run("impossible_constraints", test_impossible_constraints);
    run("cancellation_and_budget", test_cancellation_and_budget);
    run("vertical_faces", test_vertical_faces);
    run("validate_detects_tampering", test_validate_detects_tampering);
    run("crossings", test_crossings);
    run("review_regressions", test_review_regressions);
    std::printf("%d checks, %d failures\n", checks, failures);
    return failures == 0 ? 0 : 1;
}
