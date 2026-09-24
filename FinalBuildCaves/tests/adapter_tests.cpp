// MOOCoW sectioned-CSV adapter tests (fbs::caves::to_moocow_csv).
//
// Normal build: linked against fbs::caves (real validate/stitch/generate).
// Standalone build: -DFBS_ADAPTER_STUB_REGION with tests/adapter_region_stub.cpp.
//
//   caves_adapter_test                        run all native checks
//   caves_adapter_test --write-fixtures DIR   also write CSV fixtures for the
//                                             unmodified-MOOCoW .NET acceptance
//                                             (tests/adapter_acceptance.py)

#include "fbs/caves.hpp"

#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <sstream>
#include <string>
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
            std::cerr << __FILE__ << ":" << __LINE__ << ": CHECK failed: " #cond "\n";              \
        }                                                                                            \
    } while (0)

#define CHECK_EQ(a, b)                                                                               \
    do {                                                                                             \
        ++checks;                                                                                    \
        const auto& va_ = (a);                                                                       \
        const auto& vb_ = (b);                                                                       \
        if (!(va_ == vb_)) {                                                                         \
            ++failures;                                                                              \
            std::cerr << __FILE__ << ":" << __LINE__ << ": CHECK_EQ failed: " #a " == " #b "\n  left:  " \
                      << va_ << "\n  right: " << vb_ << "\n";                                        \
        }                                                                                            \
    } while (0)

const std::string map_id = "5f0c7d2e-9a41-4b8e-8d3c-2e7f61a0b9c4";

// ---------------------------------------------------------------- reader
// Mirrors CsvImport.ParseSections / ParseCsvLine so tests read the output the
// same way MOOCoW does (line split first, trimmed lines, quote toggling).

using Row = std::vector<std::string>;
using Sections = std::map<std::string, std::vector<Row>>;

Row parse_line(const std::string& line) {
    Row out;
    std::string cur;
    bool quoted = false;
    for (std::size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];
        if (c == '"') {
            if (quoted && i + 1 < line.size() && line[i + 1] == '"') {
                cur += '"';
                ++i;
            } else {
                quoted = !quoted;
            }
        } else if (c == ',' && !quoted) {
            out.push_back(cur);
            cur.clear();
        } else {
            cur += c;
        }
    }
    if (quoted) throw std::runtime_error("unclosed quote");
    out.push_back(cur);
    return out;
}

std::string trim(const std::string& s) {
    const auto b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return {};
    const auto e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}

Sections parse(const std::string& text, std::vector<std::string>* order = nullptr) {
    Sections sections;
    std::string current;
    std::istringstream in(text);
    std::string raw;
    while (std::getline(in, raw)) {
        const std::string line = trim(raw);
        if (line.empty()) continue;
        if (line.size() >= 6 && line.compare(0, 3, "===") == 0 && line.compare(line.size() - 3, 3, "===") == 0) {
            current = trim(line.substr(3, line.size() - 6));
            if (sections.count(current)) throw std::runtime_error("duplicate section");
            sections[current];
            if (order) order->push_back(current);
        } else if (!current.empty()) {
            sections[current].push_back(parse_line(line));
        }
    }
    return sections;
}

const Row* find_row(const Sections& s, const std::string& section, const std::string& id) {
    const auto it = s.find(section);
    if (it == s.end()) return nullptr;
    for (std::size_t i = 1; i < it->second.size(); ++i)
        if (!it->second[i].empty() && it->second[i][0] == id) return &it->second[i];
    return nullptr;
}

std::vector<const Row*> child_rows(const Sections& s, const std::string& section, const std::string& owner) {
    std::vector<const Row*> out;
    const auto it = s.find(section);
    if (it == s.end()) return out;
    for (std::size_t i = 1; i < it->second.size(); ++i)
        if (it->second[i][0] == owner) out.push_back(&it->second[i]);
    return out;
}

double num(const std::string& s) {
    char* end = nullptr;
    const double v = std::strtod(s.c_str(), &end);
    if (end == s.c_str() || *end != '\0') throw std::runtime_error("bad number '" + s + "'");
    return v;
}

bool same_bits(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }

// ---------------------------------------------------------------- fixture
// Hand-constructed region: irregular concave outline, sloped floor, shelf
// with a climbable cliff and a void cliff, required/optional/sealed ports,
// quoted and comma-bearing text, JSON-hostile parameters and tags, extreme
// int32 seeds, locked geometry and a variable tunnel profile. All local
// coordinates are dyadic so world = origin + local is exact both ways.

const Vec3 region_size{64, 64, 32};

Region constructed(RegionKey key = {-3, 2, -1}) {
    Region r;
    r.id = "8d8f1b8e-2a6c-4c55-9b0f-3f9e4d1a7c01";
    r.settings_id = "6e1f0a52-7c3d-4b9e-a2f4-0d8c5b7e3a91";
    r.world_seed = 0xfeedfacecafebeefull;
    r.key = key;
    r.size = region_size;
    r.origin = {static_cast<double>(key.x) * region_size.x, static_cast<double>(key.y) * region_size.y,
                static_cast<double>(key.z) * region_size.z};

    Cavern hall;
    hall.id = "1a2b3c4d-0000-4000-8000-00000000c001";
    hall.code = "C-1";
    hall.name = "Hall, \"Grand\"";
    hall.center = {16, 15, 6};
    hall.bounds = {{2, 2.5, 1}, {30, 28, 20}};
    hall.boundary = {{4, 4}, {28, 2.5}, {30, 14}, {20, 12}, {22, 26}, {6, 28}, {2, 16}};
    hall.authoring = Authoring::exact;
    hall.locked = true;
    hall.reserved_clearance = 1.5;
    hall.generation.seed = std::numeric_limits<std::int32_t>::min();
    hall.generation.parameters = {{"note", "line1\nline2"}, {"style", "wet, \"dripping\""},
                                  {"unicode", "\xC3\xA9\xF0\x9F\x98\x80<&>'+`\\"}, {"empty", ""}};
    hall.tags = {"wet", "Authored", "has,comma", "quote\"tag"};

    Floor lower;
    lower.id = "1a2b3c4d-0000-4000-8000-00000000f001";
    lower.cavern_id = hall.id;
    lower.code = "F-1";
    lower.name = "Lower floor";
    lower.boundary = {{5, 5}, {18, 4}, {20, 11}, {8, 14}};
    lower.base_z = 2;
    lower.slope_x = 0.125;
    lower.slope_y = -0.0625;
    lower.variation_amplitude = 0.25;
    lower.variation_seed = -7;
    lower.material = "lime, stone";
    lower.traversal = Traversal::walkable;
    lower.tags = {"main"};

    Floor shelf;
    shelf.id = "1a2b3c4d-0000-4000-8000-00000000f002";
    shelf.cavern_id = hall.id;
    shelf.code = "F-2";
    shelf.name = "Upper shelf";
    shelf.boundary = {{20, 14}, {27, 15}, {26, 24}, {19, 23}};
    shelf.base_z = 6;
    shelf.material = "limestone";
    shelf.traversal = Traversal::difficult;

    Cliff climb;
    climb.id = "1a2b3c4d-0000-4000-8000-00000000e001";
    climb.cavern_id = hall.id;
    climb.from_floor_id = lower.id;
    climb.to_floor_id = shelf.id;
    climb.code = "cliff-1";
    climb.boundary = {{18, 12}, {20, 13.5}, {21, 14}};
    climb.transition = Transition::climbable_cliff;
    climb.height = 4;

    Cliff drop;
    drop.id = "1a2b3c4d-0000-4000-8000-00000000e002";
    drop.cavern_id = hall.id;
    drop.from_floor_id = shelf.id;
    drop.code = "cliff-2";
    drop.boundary = {{27, 15}, {26, 24}};
    drop.transition = Transition::impassable_cliff;
    drop.height = 6;

    Port gate;
    gate.id = "1a2b3c4d-0000-4000-8000-00000000a001";
    gate.cavern_id = hall.id;
    gate.floor_id = lower.id;
    gate.code = "P-1";
    gate.name = "South gate";
    // On the sloped floor: 2 + 0.125*(10-5) - 0.0625*(6-5).
    gate.position = {10, 6, 2.5625};
    gate.facing = {0, -1, 0};
    gate.width = 4;
    gate.height = 3;
    gate.state = PortState::required;
    gate.traversal = Traversal::walkable;

    Cavern pocket;
    pocket.id = "1a2b3c4d-0000-4000-8000-00000000c002";
    pocket.code = "C-2";
    pocket.name = "Side pocket";
    pocket.center = {50, 50, 6};
    pocket.bounds = {{40, 40, 2}, {60, 60, 14}};
    pocket.authoring = Authoring::constrained_procedural;
    pocket.generation.seed = 42;
    pocket.tags = {"pocket"};

    Floor pocket_floor;
    pocket_floor.id = "1a2b3c4d-0000-4000-8000-00000000f003";
    pocket_floor.cavern_id = pocket.id;
    pocket_floor.code = "F-3";
    pocket_floor.name = "Pocket floor";
    pocket_floor.boundary = {{42, 42}, {58, 42}, {58, 58}, {42, 58}};
    pocket_floor.base_z = 4;

    Port side;
    side.id = "1a2b3c4d-0000-4000-8000-00000000a002";
    side.cavern_id = pocket.id;
    side.floor_id = pocket_floor.id;
    side.code = "P-2";
    side.name = "West mouth";
    side.position = {45, 50, 4};
    side.facing = {-1, 0, 0};
    side.width = 5;
    side.height = 4;
    side.state = PortState::optional;

    Port sealed;
    sealed.id = "1a2b3c4d-0000-4000-8000-00000000a003";
    sealed.cavern_id = pocket.id;
    sealed.code = "P-3";
    sealed.name = "Sealed vent";
    sealed.position = {60, 55, 8};
    sealed.facing = {1, 0, 0};
    sealed.width = 2;
    sealed.height = 2;
    sealed.state = PortState::sealed;
    sealed.traversal = Traversal::impassable;

    Tunnel crawl;
    crawl.id = "1a2b3c4d-0000-4000-8000-00000000d001";
    crawl.start_port_id = gate.id;
    crawl.end_port_id = side.id;
    crawl.code = "T-1";
    crawl.name = "Crawl, the \"long\" one";
    crawl.centerline = {{gate.position, 4, 3}, {{34, 20, 3}, 5.5, 3.75}, {side.position, 5, 4}};
    crawl.maximum_slope_degrees = 30;
    crawl.minimum_clearance = 3;
    crawl.traversal = Traversal::hazardous;
    crawl.locked = true;
    crawl.generation.seed = std::numeric_limits<std::int32_t>::max();
    crawl.generation.parameters = {{"route", "hand, \"made\""}};
    crawl.tags = {"loop", "Authored"};
    crawl.source_node_ids = {7, 9, 11};

    r.caverns = {hall, pocket};
    r.floors = {lower, shelf, pocket_floor};
    r.cliffs = {climb, drop};
    r.ports = {gate, side, sealed};
    r.tunnels = {crawl};
    return r;
}

// A small isolated neighbour with distinct IDs and codes.
Region neighbour(RegionKey key) {
    Region r;
    r.id = "8d8f1b8e-2a6c-4c55-9b0f-3f9e4d1a7c02";
    r.settings_id = "6e1f0a52-7c3d-4b9e-a2f4-0d8c5b7e3a91";
    r.world_seed = 0xfeedfacecafebeefull;
    r.key = key;
    r.size = region_size;
    r.origin = {static_cast<double>(key.x) * region_size.x, static_cast<double>(key.y) * region_size.y,
                static_cast<double>(key.z) * region_size.z};
    Cavern c;
    c.id = "2b3c4d5e-0000-4000-8000-00000000c003";
    c.code = "C-3";
    c.name = "Neighbour vault";
    c.center = {32, 32, 8};
    c.bounds = {{24, 24, 4}, {40, 40, 16}};
    c.boundary = {{24, 24}, {40, 26}, {38, 40}, {26, 38}};
    c.tags = {"isolated"};
    // Case-sensitive neighbour of MOOCoW's legacy "parameters" wrapper key.
    c.generation.parameters = {{"Parameters", "kept"}, {"z", "\xE2\x82\xAC"}, {"\xEF\xBC\xA1", "fullwidth"},
                               {"\xF0\x9F\x98\x80", "astral"}};
    Floor f;
    f.id = "2b3c4d5e-0000-4000-8000-00000000f004";
    f.cavern_id = c.id;
    f.code = "F-4";
    f.name = "Vault floor";
    f.boundary = {{26, 26}, {38, 27}, {36, 38}, {27, 36}};
    f.base_z = 5;
    r.caverns = {c};
    r.floors = {f};
    return r;
}

template <class F>
bool rejects(F&& f, ErrorCode* code = nullptr, std::string* message = nullptr) {
    try {
        f();
    } catch (const Error& e) {
        if (code) *code = e.code;
        if (message) *message = e.what();
        return true;
    }
    return false;
}

#define CHECK_REJECTS(expr)                                                                          \
    do {                                                                                             \
        std::string message_;                                                                        \
        const bool rejected_ = rejects([&] { (void)(expr); }, nullptr, &message_);                   \
        ++checks;                                                                                    \
        if (!rejected_) {                                                                            \
            ++failures;                                                                              \
            std::cerr << __FILE__ << ":" << __LINE__ << ": expected rejection: " #expr "\n";         \
        }                                                                                            \
    } while (0)

// With the stub, the adapter's own message is asserted. With the real
// library, validate() may legitimately reject first in its own words, so only
// the rejection itself is asserted (the stub build covers the adapter texts).
#ifdef FBS_ADAPTER_STUB_REGION
constexpr bool exact_messages = true;
#else
constexpr bool exact_messages = false;
#endif

#define CHECK_REJECTS_WITH(expr, needle)                                                             \
    do {                                                                                             \
        std::string message_;                                                                        \
        const bool rejected_ = rejects([&] { (void)(expr); }, nullptr, &message_);                   \
        ++checks;                                                                                    \
        if (!rejected_ || (exact_messages && message_.find(needle) == std::string::npos)) {          \
            ++failures;                                                                              \
            std::cerr << __FILE__ << ":" << __LINE__ << ": expected rejection mentioning '" << needle \
                      << "': " #expr "\n  got: " << (rejected_ ? message_ : "(accepted)") << "\n";   \
        }                                                                                            \
    } while (0)

// ---------------------------------------------------------------- tests

const char* expected_headers[][2] = {
    {"MAP", "id,code,name,min_x,min_y,min_z,max_x,max_y,max_z"},
    {"ELEVATION_BANDS", "id,map_id,code,name,min_z,max_z,display_order"},
    {"CAVERNS", "id,map_id,code,name,center_x,center_y,center_z,min_x,min_y,min_z,max_x,max_y,max_z,"
                "authoring_mode,is_locked,reserved_clearance_m,generation_seed,tags,generation_json"},
    {"CAVERN_OUTLINE_POINTS", "cavern_id,point_index,x,y"},
    {"FLOOR_REGIONS", "id,cavern_id,code,name,base_z,slope_x,slope_y,variation_amplitude_m,variation_seed,"
                      "material_code,traversal_class,tags"},
    {"FLOOR_REGION_POINTS", "floor_region_id,point_index,x,y"},
    {"CLIFF_EDGES", "id,cavern_id,from_floor_region_id,to_floor_region_id,code,transition_kind,height_m"},
    {"CLIFF_EDGE_POINTS", "cliff_edge_id,point_index,x,y"},
    {"CAVERN_PORTS", "id,cavern_id,floor_region_id,code,name,position_x,position_y,position_z,"
                     "facing_x,facing_y,facing_z,width_m,height_m,port_state,traversal_class"},
    {"TUNNELS", "id,map_id,start_port_id,end_port_id,code,name,maximum_slope_degrees,minimum_clearance_m,"
                "traversal_class,is_locked,generation_seed,tags,generation_json"},
    {"TUNNEL_POINTS", "tunnel_id,point_index,center_x,center_y,center_z,width_m,height_m"},
};

std::string join(const Row& r) {
    std::string s;
    for (std::size_t i = 0; i < r.size(); ++i) s += (i ? "," : "") + r[i];
    return s;
}

void test_layout_and_values() {
    const Region r = constructed();
    const std::string text = to_moocow_csv({r}, map_id, "Adapter, \"constructed\" map");
    std::vector<std::string> order;
    const Sections s = parse(text, &order);

    // Exact section order and headers of CsvExport; no extra sections.
    CHECK_EQ(order.size(), std::size(expected_headers));
    for (std::size_t i = 0; i < std::size(expected_headers) && i < order.size(); ++i) {
        CHECK_EQ(order[i], std::string(expected_headers[i][0]));
        CHECK_EQ(join(s.at(order[i]).at(0)), std::string(expected_headers[i][1]));
    }
    // CsvExport line structure: LF only, one blank line between sections.
    CHECK(text.find('\r') == std::string::npos);
    CHECK(text.compare(0, 14, "=== MAP ===\nid") == 0);
    CHECK(text.find("\n\n=== ELEVATION_BANDS ===\n") != std::string::npos);
    CHECK(text.back() == '\n' && text[text.size() - 2] != '\n');

    const Vec3 o = r.origin;
    CHECK_EQ(o.x, -192.0);
    CHECK_EQ(o.z, -32.0);

    // MAP: fresh derived code, quoted name, bounds contain the whole region box.
    const Row& m = s.at("MAP").at(1);
    CHECK_EQ(m.size(), std::size_t(9));
    CHECK_EQ(m[0], map_id);
    CHECK_EQ(m[1], "fbs-caves-" + map_id);
    CHECK_EQ(m[2], std::string("Adapter, \"constructed\" map"));
    CHECK(num(m[3]) <= o.x && num(m[4]) <= o.y && num(m[5]) <= o.z);
    CHECK(num(m[6]) >= o.x + 64 && num(m[7]) >= o.y + 64 && num(m[8]) >= o.z + 32);
    CHECK(text.find(",\"Adapter, \"\"constructed\"\" map\",") != std::string::npos);

    // One band for the single Z layer, covering the map Z range.
    CHECK_EQ(s.at("ELEVATION_BANDS").size(), std::size_t(2));
    const Row& band = s.at("ELEVATION_BANDS").at(1);
    CHECK_EQ(band[1], map_id);
    CHECK_EQ(band[2], std::string("fbs-z-1"));
    CHECK_EQ(num(band[4]), num(m[5]));
    CHECK_EQ(num(band[5]), num(m[8]));
    CHECK_EQ(band[6], std::string("0"));

    // Cavern row: world offsets, exact authored metadata, canonical tags/JSON.
    const Row* hall = find_row(s, "CAVERNS", r.caverns[0].id);
    CHECK(hall != nullptr);
    if (hall) {
        CHECK_EQ(hall->size(), std::size_t(19));
        CHECK_EQ((*hall)[1], map_id);
        CHECK_EQ((*hall)[2], std::string("C-1"));
        CHECK_EQ((*hall)[3], std::string("Hall, \"Grand\""));
        CHECK_EQ(num((*hall)[4]) - o.x, 16.0);
        CHECK_EQ(num((*hall)[5]) - o.y, 15.0);
        CHECK_EQ(num((*hall)[6]) - o.z, 6.0);
        CHECK_EQ(num((*hall)[8]) - o.y, 2.5);
        CHECK_EQ((*hall)[13], std::string("exact"));
        CHECK_EQ((*hall)[14], std::string("1"));
        CHECK_EQ((*hall)[15], std::string("1.5"));
        CHECK_EQ((*hall)[16], std::string("-2147483648"));
        CHECK_EQ((*hall)[17], std::string("[\"Authored\",\"has,comma\",\"quote\\u0022tag\",\"wet\"]"));
        CHECK_EQ((*hall)[18], std::string("{\"empty\":\"\",\"note\":\"line1\\nline2\",\"style\":\"wet, "
                                          "\\u0022dripping\\u0022\",\"unicode\":\"\\u00E9\\uD83D\\uDE00"
                                          "\\u003C\\u0026\\u003E\\u0027\\u002B\\u0060\\\\\"}"));
    }
    const std::string hall_line =
        "1a2b3c4d-0000-4000-8000-00000000c001," + map_id +
        ",C-1,\"Hall, \"\"Grand\"\"\",-176,143,-26,-190,130.5,-31,-162,156,-12,exact,1,1.5,-2147483648,"
        "\"[\"\"Authored\"\",\"\"has,comma\"\",\"\"quote\\u0022tag\"\",\"\"wet\"\"]\","
        "\"{\"\"empty\"\":\"\"\"\",\"\"note\"\":\"\"line1\\nline2\"\",\"\"style\"\":\"\"wet, \\u0022dripping\\u0022\"\","
        "\"\"unicode\"\":\"\"\\u00E9\\uD83D\\uDE00\\u003C\\u0026\\u003E\\u0027\\u002B\\u0060\\\\\"\"}\"\n";
    CHECK(text.find(hall_line) != std::string::npos);

    const Row* pocket = find_row(s, "CAVERNS", r.caverns[1].id);
    CHECK(pocket != nullptr);
    if (pocket) {
        CHECK_EQ((*pocket)[13], std::string("constrained_procedural"));
        CHECK_EQ((*pocket)[14], std::string("0"));
        CHECK_EQ((*pocket)[17], std::string("[\"pocket\"]"));
        CHECK_EQ((*pocket)[18], std::string("{}"));
    }

    // Outline: every irregular vertex, contiguous indexes, world XY; none for pocket.
    const auto outline = child_rows(s, "CAVERN_OUTLINE_POINTS", r.caverns[0].id);
    CHECK_EQ(outline.size(), r.caverns[0].boundary.size());
    for (std::size_t i = 0; i < outline.size(); ++i) {
        CHECK_EQ((*outline[i])[1], std::to_string(i));
        CHECK(same_bits(num((*outline[i])[2]) - o.x, r.caverns[0].boundary[i].x));
        CHECK(same_bits(num((*outline[i])[3]) - o.y, r.caverns[0].boundary[i].y));
    }
    CHECK(child_rows(s, "CAVERN_OUTLINE_POINTS", r.caverns[1].id).empty());

    // Sloped floor: gradients unchanged, base_z offset, quoted material, seed.
    const Row* lower = find_row(s, "FLOOR_REGIONS", r.floors[0].id);
    CHECK(lower != nullptr);
    if (lower) {
        CHECK_EQ(join(*lower), r.floors[0].id + "," + r.caverns[0].id +
                                   ",F-1,Lower floor,-30,0.125,-0.0625,0.25,-7,lime, stone,walkable,[\"main\"]");
    }
    CHECK(text.find(",-30,0.125,-0.0625,0.25,-7,\"lime, stone\",walkable,\"[\"\"main\"\"]\"\n") != std::string::npos);
    const Row* shelf = find_row(s, "FLOOR_REGIONS", r.floors[1].id);
    CHECK(shelf && (*shelf)[10] == "difficult" && (*shelf)[11] == "[]");
    CHECK_EQ(child_rows(s, "FLOOR_REGION_POINTS", r.floors[0].id).size(), std::size_t(4));

    // Cliffs: destination and void edge (empty to_floor), transition, height, path.
    const Row* climb = find_row(s, "CLIFF_EDGES", r.cliffs[0].id);
    CHECK(climb && (*climb)[3] == r.floors[1].id && (*climb)[5] == "climbable_cliff" && (*climb)[6] == "4");
    const Row* drop = find_row(s, "CLIFF_EDGES", r.cliffs[1].id);
    CHECK(drop && (*drop)[3].empty() && (*drop)[5] == "impassable_cliff" && (*drop)[6] == "6");
    const auto climb_points = child_rows(s, "CLIFF_EDGE_POINTS", r.cliffs[0].id);
    CHECK_EQ(climb_points.size(), std::size_t(3));
    if (climb_points.size() == 3) CHECK_EQ(num((*climb_points[1])[3]) - o.y, 13.5);

    // Ports: states, optional floor, facing, dimensions, world position.
    const Row* gate = find_row(s, "CAVERN_PORTS", r.ports[0].id);
    CHECK(gate != nullptr);
    if (gate) {
        CHECK_EQ(join(*gate), r.ports[0].id + "," + r.caverns[0].id + "," + r.floors[0].id +
                                  ",P-1,South gate,-182,134,-29.4375,0,-1,0,4,3,required,walkable");
    }
    const Row* sealed = find_row(s, "CAVERN_PORTS", r.ports[2].id);
    CHECK(sealed && (*sealed)[2].empty() && (*sealed)[13] == "sealed" && (*sealed)[14] == "impassable");

    // Tunnel: lock, extreme seed, tags, params and variable profile.
    const Row* crawl = find_row(s, "TUNNELS", r.tunnels[0].id);
    CHECK(crawl != nullptr);
    if (crawl) {
        CHECK_EQ((*crawl)[1], map_id);
        CHECK_EQ((*crawl)[5], std::string("Crawl, the \"long\" one"));
        CHECK_EQ((*crawl)[6], std::string("30"));
        CHECK_EQ((*crawl)[7], std::string("3"));
        CHECK_EQ((*crawl)[8], std::string("hazardous"));
        CHECK_EQ((*crawl)[9], std::string("1"));
        CHECK_EQ((*crawl)[10], std::string("2147483647"));
        CHECK_EQ((*crawl)[11], std::string("[\"Authored\",\"loop\"]"));
        CHECK_EQ((*crawl)[12], std::string("{\"route\":\"hand, \\u0022made\\u0022\"}"));
    }
    const auto profile = child_rows(s, "TUNNEL_POINTS", r.tunnels[0].id);
    CHECK_EQ(profile.size(), std::size_t(3));
    for (std::size_t i = 0; i < profile.size(); ++i) {
        const Profile& p = r.tunnels[0].centerline[i];
        CHECK(same_bits(num((*profile[i])[2]) - o.x, p.center.x));
        CHECK(same_bits(num((*profile[i])[3]) - o.y, p.center.y));
        CHECK(same_bits(num((*profile[i])[4]) - o.z, p.center.z));
        CHECK(same_bits(num((*profile[i])[5]), p.width));
        CHECK(same_bits(num((*profile[i])[6]), p.height));
    }
    // Records are in MOOCoW's storage order (ORDER BY id) so the CSV equals
    // MOOCoW's own export of the stored map.
    for (const char* section : {"CAVERNS", "TUNNELS"}) {
        const auto& rows = s.at(section);
        for (std::size_t i = 2; i < rows.size(); ++i) CHECK(rows[i - 1][0] < rows[i][0]);
    }
    for (const char* section : {"FLOOR_REGIONS", "CLIFF_EDGES", "CAVERN_PORTS"}) {
        const auto& rows = s.at(section);
        for (std::size_t i = 2; i < rows.size(); ++i)
            CHECK(rows[i - 1][1] < rows[i][1] || (rows[i - 1][1] == rows[i][1] && rows[i - 1][0] < rows[i][0]));
    }
    // Provenance never becomes a column or a value.
    CHECK(text.find("source_node") == std::string::npos);
    CHECK(text.find("6e1f0a52-7c3d-4b9e-a2f4-0d8c5b7e3a91") == std::string::npos);
}

void test_determinism_and_order() {
    const Region a = constructed({-3, 2, -1});
    const Region b = neighbour({-2, 2, -1});
    const std::string ab = to_moocow_csv({a, b}, map_id);
    const std::string ba = to_moocow_csv({b, a}, map_id);
    CHECK_EQ(ab, ba);
    CHECK_EQ(ab, to_moocow_csv({a, b}, map_id));
    const Sections s = parse(ab);
    CHECK_EQ(s.at("CAVERNS").size(), std::size_t(4));
    CHECK_EQ(s.at("MAP").at(1)[2], std::string("Generated caves"));
    // Neighbour offset by its own origin.
    const Row* vault = find_row(s, "CAVERNS", b.caverns[0].id);
    CHECK(vault && num((*vault)[4]) == b.origin.x + 32 && num((*vault)[6]) == b.origin.z + 8);
    // Map spans both boxes.
    const Row& m = s.at("MAP").at(1);
    CHECK(num(m[3]) <= a.origin.x && num(m[6]) >= b.origin.x + 64);
    // Uppercase map IDs are canonicalised without changing the UUID value.
    std::string upper = map_id;
    for (auto& c : upper) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    CHECK_EQ(to_moocow_csv({a, b}, upper), ab);
}

void test_multi_layer_bands() {
    const Region low = constructed({0, 0, 0});
    const Region high = neighbour({0, 0, 1});
    const Sections s = parse(to_moocow_csv({high, low}, map_id));
    const auto& bands = s.at("ELEVATION_BANDS");
    CHECK_EQ(bands.size(), std::size_t(3));
    if (bands.size() == 3) {
        CHECK_EQ(bands[1][2], std::string("fbs-z0"));
        CHECK_EQ(bands[2][2], std::string("fbs-z1"));
        CHECK_EQ(num(bands[1][5]), 32.0);  // touching, never exclusively overlapping
        CHECK_EQ(num(bands[2][4]), 32.0);
        CHECK_EQ(bands[1][6], std::string("0"));
        CHECK_EQ(bands[2][6], std::string("1"));
        CHECK(bands[1][0] != bands[2][0]);
    }
}

void test_precision_contract() {
    // Non-dyadic local values with a non-integer origin: world = fl(o + l),
    // so world - o recovers l to within half an ulp of the world value.
    Region r = constructed({0, 0, 0});
    r.origin = {1000.1, -2000.3, 0.7};
    r.caverns[1].center = {50.1, 49.9, 6.3};
    const Sections s = parse(to_moocow_csv({r}, map_id));
    const Row* pocket = find_row(s, "CAVERNS", r.caverns[1].id);
    CHECK(pocket != nullptr);
    if (pocket) {
        const double local[3] = {50.1, 49.9, 6.3};
        const double origin[3] = {r.origin.x, r.origin.y, r.origin.z};
        for (int a = 0; a < 3; ++a) {
            const double world = num((*pocket)[4 + a]);
            CHECK(same_bits(world, origin[a] + local[a]));  // G17 round-trips the double exactly
            const double half_ulp = (std::nextafter(std::fabs(world), INFINITY) - std::fabs(world)) / 2;
            CHECK(std::fabs((world - origin[a]) - local[a]) <= half_ulp);
        }
    }
}

void test_parameters_key_rejected() {
    // MOOCoW cannot store and reload a generation parameter named "parameters"
    // (see adapter_acceptance.cs probe); the adapter refuses it explicitly.
    Region r = constructed();
    r.tunnels[0].generation.parameters["parameters"] = "literal";
    CHECK_REJECTS_WITH(to_moocow_csv({r}, map_id), "named 'parameters'");
    r = constructed();
    r.caverns[0].generation.parameters["parameters"] = "{}";
    CHECK_REJECTS_WITH(to_moocow_csv({r}, map_id), "named 'parameters'");
    // Only the exact, case-sensitive key is special in MOOCoW.
    r = constructed();
    r.caverns[0].generation.parameters["Parameters"] = "kept";
    r.caverns[0].generation.parameters["parameters.v2"] = "kept";
    const Sections s = parse(to_moocow_csv({r}, map_id));
    const Row* hall = find_row(s, "CAVERNS", r.caverns[0].id);
    CHECK(hall && (*hall)[18].find("\"Parameters\":\"kept\"") != std::string::npos &&
          (*hall)[18].find("\"parameters.v2\":\"kept\"") != std::string::npos);
}

void test_rejections() {
    const Region good = constructed();
    ErrorCode code{};

    // Map identity and input shape.
    CHECK_REJECTS_WITH(to_moocow_csv({good}, ""), "map_id");
    CHECK_REJECTS_WITH(to_moocow_csv({good}, "not-a-uuid"), "map_id");
    CHECK_REJECTS_WITH(to_moocow_csv({good}, "00000000-0000-0000-0000-000000000000"), "nil");
    CHECK_REJECTS(to_moocow_csv({good}, "{5f0c7d2e-9a41-4b8e-8d3c-2e7f61a0b9c4}"));
    CHECK_REJECTS(to_moocow_csv({good}, "5f0c7d2e9a414b8e8d3c2e7f61a0b9c4"));
    CHECK_REJECTS(to_moocow_csv({good}, "5f0c7d2e-9a41-4b8e-8d3c-2e7f61a0b9cg"));
    CHECK_REJECTS_WITH(to_moocow_csv({good}, map_id, ""), "map_name");
    CHECK_REJECTS_WITH(to_moocow_csv({good}, map_id, " padded"), "whitespace");
    CHECK_REJECTS_WITH(to_moocow_csv({good}, map_id, "two\nlines"), "line break");
    CHECK(rejects([&] { to_moocow_csv({}, map_id); }, &code) && code == ErrorCode::invalid_input);
    Region empty = good;
    empty.caverns.clear();
    empty.floors.clear();
    empty.cliffs.clear();
    empty.ports.clear();
    empty.tunnels.clear();
    CHECK_REJECTS_WITH(to_moocow_csv({empty}, map_id), "no caverns");
    CHECK_REJECTS_WITH(to_moocow_csv({good}, good.caverns[0].id), "duplicate entity ID");

    // Duplicate records: never silently deduplicated.
    CHECK_REJECTS_WITH(to_moocow_csv({good, good}, map_id), "more than once");
    Region twin = neighbour({-2, 2, -1});
    twin.caverns[0].id = good.caverns[0].id;
    twin.floors[0].cavern_id = good.caverns[0].id;
    CHECK_REJECTS_WITH(to_moocow_csv({good, twin}, map_id), "duplicate entity ID");
    Region same_floor = good;
    same_floor.floors[1].id = same_floor.floors[0].id;
    CHECK_REJECTS_WITH(to_moocow_csv({same_floor}, map_id), "duplicate entity ID");
    Region upper_dup = good;
    std::string upper_id = good.ports[0].id;
    for (auto& c : upper_id) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    upper_dup.ports[2].id = upper_id;
    CHECK_REJECTS_WITH(to_moocow_csv({upper_dup}, map_id), "duplicate entity ID");
    Region code_clash = neighbour({-2, 2, -1});
    code_clash.caverns[0].code = "c-1";  // map-scope, case-insensitive
    CHECK_REJECTS_WITH(to_moocow_csv({good, code_clash}, map_id), "duplicated in map:");
    Region band_clash = good;
    band_clash.caverns[1].code = "FBS-Z-1";
    CHECK_REJECTS_WITH(to_moocow_csv({band_clash}, map_id), "duplicated in map:");
    Region port_code = good;
    port_code.ports[0].code = "f-1";  // same cavern scope as floor F-1
    CHECK_REJECTS_WITH(to_moocow_csv({port_code}, map_id), "duplicated in cavern:");

    // Inconsistent region settings and lattice.
    Region other_seed = neighbour({-2, 2, -1});
    other_seed.world_seed ^= 1;
    CHECK(rejects([&] { to_moocow_csv({good, other_seed}, map_id); }, &code) &&
          code == ErrorCode::incompatible_boundary);
    Region other_settings = neighbour({-2, 2, -1});
    other_settings.settings_id = "6e1f0a52-7c3d-4b9e-a2f4-0d8c5b7e3a92";
    CHECK_REJECTS_WITH(to_moocow_csv({good, other_settings}, map_id), "settings_id");
    Region other_size = neighbour({-2, 2, -1});
    other_size.size.x = 65;
    CHECK_REJECTS_WITH(to_moocow_csv({good, other_size}, map_id), "size");
    Region off_lattice = neighbour({-2, 2, -1});
    off_lattice.origin.x += 0.5;
    CHECK_REJECTS_WITH(to_moocow_csv({good, off_lattice}, map_id), "lattice");
    Region old = good;
    old.version = schema_version + 1;
    CHECK(rejects([&] { to_moocow_csv({old}, map_id); }, &code) && code == ErrorCode::incompatible_version);
    Region other_algorithm = good;
    other_algorithm.algorithm = "something-else";
    CHECK(rejects([&] { to_moocow_csv({other_algorithm}, map_id); }, &code) &&
          code == ErrorCode::incompatible_version);
    Region bad_size = good;
    bad_size.size.z = 0;
    CHECK_REJECTS_WITH(to_moocow_csv({bad_size}, map_id), "size");

    // Non-finite transport numbers anywhere.
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    Region v = good;
    v.tunnels[0].centerline[1].center.z = nan;
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "finite");
    v = good;
    v.origin.y = inf;
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "finite");
    v = good;
    v.floors[0].slope_x = nan;
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "finite");
    v = good;
    v.ports[1].facing.y = inf;
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "finite");
    v = good;
    v.caverns[0].reserved_clearance = -inf;
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "finite");
    v = good;
    v.origin.x = 1.7e308;
    v.caverns[0].center.x = 1.7e308;
    CHECK_REJECTS(to_moocow_csv({v}, map_id));

    // Values the existing importer/domain would change or drop.
    v = good;
    v.caverns[0].name = "Hall\r\nTwo";
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "line break");
    v = good;
    v.floors[0].material = "limestone ";
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "whitespace");
    v = good;
    v.ports[0].code = "\xC2\xA0P-1";  // NBSP is .NET whitespace
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "whitespace");
    v = good;
    v.caverns[0].tags.push_back("WET");
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "ignoring case");
    v = good;
    v.caverns[0].tags.push_back("h\xC3\xB6hle");
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "non-ASCII");
    v = good;
    v.tunnels[0].tags.push_back("");
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "empty");
    v = good;
    v.tunnels[0].generation.parameters[" padded"] = "x";
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "whitespace");
    v = good;
    v.caverns[0].name = "bad \xC3\x28 utf8";
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "UTF-8");
    v = good;
    v.ports[1].facing = {-2, 0, 0};
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "normalized");
    v = good;
    v.ports[1].facing = {0, 0, 0};
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "zero magnitude");
    v = good;
    v.tunnels[0].id = "not-a-uuid";
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "UUID");
    v = good;
    v.ports[0].floor_id = "00000000-0000-0000-0000-000000000000";
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "nil");

    // MOOCoW validator rules.
    v = good;
    v.tunnels.clear();
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "required but no tunnel");
    v = good;
    v.tunnels[0].end_port_id = v.ports[2].id;
    v.tunnels[0].centerline.back().center = v.ports[2].position;
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "sealed");
    v = good;
    v.tunnels[0].maximum_slope_degrees = 0.5;
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "slope");
    v = good;
    v.tunnels[0].maximum_slope_degrees = 90;
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "[0, 90)");
    v = good;
    v.tunnels[0].centerline[1].height = 2.5;
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "clearance");
    v = good;
    v.tunnels[0].centerline.front().center.x += 10;
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "misses its start port");
    v = good;
    v.tunnels[0].end_port_id = v.tunnels[0].start_port_id;
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "two different ports");
    v = good;
    v.tunnels[0].start_port_id = "1a2b3c4d-0000-4000-8000-0000000000ff";
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "missing port");
    v = good;
    v.caverns[1].tags.clear();
    v.tunnels.clear();
    v.ports[0].state = PortState::optional;
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "not tagged isolated");
    v = good;
    v.cliffs[0].transition = Transition::unconnected;
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "isolated floor component");
    v = good;
    v.cliffs[1].transition = Transition::one_way_drop;
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "one-way drop");
    v = good;
    v.cliffs[0].height = 0;
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "positive height");
    v = good;
    v.cliffs[0].to_floor_id = v.floors[2].id;
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "another cavern");
    v = good;
    v.caverns[0].boundary = {{4, 4}, {28, 2.5}, {6, 28}, {30, 14}};  // bow-tie
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "edges intersect");
    v = good;
    v.caverns[0].boundary.clear();
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "no outline");
    v = good;
    v.caverns[0].boundary[3] = {31, 12};
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "outside its bounds");
    v = good;
    v.floors[1].boundary = {{20, 14}, {20, 14}, {26, 24}, {19, 23}};
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "consecutive duplicate");
    v = good;
    v.ports[0].position.x = 19.5;  // inside cavern bounds, outside its floor
    v.tunnels[0].centerline.front().center.x = 19.5;
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "outside its floor");
    v = good;
    v.ports[2].position.z = 20;
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "outside its cavern bounds");
    v = good;
    v.caverns[1].center.z = 20;
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "center is outside");
    v = good;
    v.floors[0].variation_amplitude = -1;
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "negative");
    v = good;
    v.ports[0].width = 0;
    CHECK_REJECTS_WITH(to_moocow_csv({v}, map_id), "positive");
}

#ifdef FBS_ADAPTER_STUB_REGION
// Seam handling with the stub stitch(): connected seams are added in world
// coordinates; one-sided or unconnected seams are rejected.
Region with_seam(Region r, Face face, Vec3 shared, const std::string& key, const std::string& port_id,
                 const std::string& cavern_id, const std::string& floor_id, Vec3 port_position) {
    Port p;
    p.id = port_id;
    p.cavern_id = cavern_id;
    p.floor_id = floor_id;
    p.code = "seam-port-" + key;
    p.name = "Seam port";
    p.position = port_position;
    p.facing = {face == Face::positive_x ? 1.0 : -1.0, 0, 0};
    p.state = PortState::required;
    r.ports.push_back(p);
    Boundary b;
    b.key = key;
    b.port_id = port_id;
    b.face = face;
    b.position = shared;
    r.boundaries.push_back(b);
    return r;
}

void test_stub_seams() {
    Region a = constructed({0, 0, 0});
    a.caverns[1].bounds.maximum.x = 64;
    a = with_seam(a, Face::positive_x, {64, 41, 6}, "seam-x0", "3c4d5e6f-0000-4000-8000-00000000a010",
                  a.caverns[1].id, a.floors[2].id, {56, 50, 4});
    Region b = neighbour({1, 0, 0});
    b.caverns[0].bounds.minimum.x = 0;
    b.caverns[0].tags.clear();  // linked through the seam, no isolation tag needed
    b = with_seam(b, Face::negative_x, {0, 41, 6}, "seam-x0", "3c4d5e6f-0000-4000-8000-00000000a011",
                  b.caverns[0].id, "", {8, 32, 8});
    const std::string text = to_moocow_csv({b, a}, map_id);
    CHECK_EQ(text, to_moocow_csv({a, b}, map_id));
    const Sections s = parse(text);
    const std::string seam_id = stable_id("stub-seam\nseam-x0");
    const Row* seam = find_row(s, "TUNNELS", seam_id);
    CHECK(seam != nullptr);
    const auto points = child_rows(s, "TUNNEL_POINTS", seam_id);
    CHECK_EQ(points.size(), std::size_t(3));
    if (points.size() == 3) {
        CHECK_EQ(num((*points[0])[2]), 56.0);   // a-local 56 + a.origin 0
        CHECK_EQ(num((*points[1])[2]), 64.0);   // shared face
        CHECK_EQ(num((*points[2])[2]), 72.0);   // b-local 8 + b.origin 64
    }

    // One-sided seam: the stub (like the real stitch) refuses it and the error propagates.
    Region lonely = neighbour({1, 0, 0});
    CHECK(rejects([&] { to_moocow_csv({a, lonely}, map_id); }));
    // A seam whose port the returned tunnel does not reach is rejected.
    Region wrong = b;
    wrong.boundaries[0].port_id = "3c4d5e6f-0000-4000-8000-00000000a0ff";
    CHECK(rejects([&] { to_moocow_csv({a, wrong}, map_id); }));
}
#else
// Real generator output, when the regional library is linked.
std::vector<Region> generated_pair() {
    Request q;
    q.world_seed = 20260924;
    q.key = {-1, 0, 0};
    Region west = generate(q);
    q.key = {0, 0, 0};
    Region east = generate(q);
    return {west, east};
}

void test_generated(std::vector<Region>* keep) {
    std::vector<Region> pair;
    try {
        pair = generated_pair();
    } catch (const std::exception& e) {
        ++failures;
        std::cerr << "generate() failed: " << e.what() << "\n";
        return;
    }
    // Equal requests reproduce equal regions (and therefore equal CSV bytes).
    try {
        const std::vector<Region> again = generated_pair();
        CHECK_EQ(to_json(again[0]), to_json(pair[0]));
        CHECK_EQ(to_json(again[1]), to_json(pair[1]));
    } catch (const std::exception& e) {
        ++failures;
        std::cerr << "second generate() failed: " << e.what() << "\n";
    }
    std::string text;
    try {
        text = to_moocow_csv(pair, map_id, "Generated pair");
    } catch (const std::exception& e) {
        ++failures;
        std::cerr << "to_moocow_csv(generated) failed: " << e.what() << "\n";
        return;
    }
    const Sections s = parse(text);
    std::size_t tunnels = 0;
    for (const auto& r : pair) tunnels += r.tunnels.size();
    CHECK(s.at("TUNNELS").size() - 1 > tunnels);  // stitched seam tunnel(s) added
    std::size_t caverns = 0;
    for (const auto& r : pair) caverns += r.caverns.size();
    CHECK_EQ(s.at("CAVERNS").size() - 1, caverns);
    const Row& m = s.at("MAP").at(1);
    for (std::size_t i = 1; i < s.at("TUNNEL_POINTS").size(); ++i) {
        const Row& p = s.at("TUNNEL_POINTS")[i];
        CHECK(num(p[2]) >= num(m[3]) && num(p[2]) <= num(m[6]));
        CHECK(num(p[4]) - num(p[6]) / 2 >= num(m[5]) && num(p[4]) + num(p[6]) / 2 <= num(m[8]));
    }
    CHECK_EQ(text, to_moocow_csv({pair[1], pair[0]}, map_id, "Generated pair"));
    if (keep) *keep = pair;
}
#endif

// MOOCoW record IDs are primary keys across the whole database, so every
// fixture that is imported next to another one needs its own record IDs.
Region with_id_group(Region r, const std::string& group) {
    auto re = [&](std::string& id) {
        if (id.size() == 36) id.replace(9, 4, group);
    };
    re(r.id);
    for (auto& c : r.caverns) re(c.id);
    for (auto& f : r.floors) { re(f.id); re(f.cavern_id); }
    for (auto& c : r.cliffs) { re(c.id); re(c.cavern_id); re(c.from_floor_id); re(c.to_floor_id); }
    for (auto& p : r.ports) { re(p.id); re(p.cavern_id); re(p.floor_id); }
    for (auto& t : r.tunnels) { re(t.id); re(t.start_port_id); re(t.end_port_id); }
    for (auto& b : r.boundaries) re(b.port_id);
    return r;
}

bool write_file(const std::string& path, const std::string& text) {
    std::ofstream out(path, std::ios::binary);
    out << text;
    return static_cast<bool>(out);
}

} // namespace

int main(int argc, char** argv) {
    std::string fixtures;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--write-fixtures" && i + 1 < argc) fixtures = argv[++i];
    }

    const std::vector<std::pair<const char*, std::function<void()>>> tests = {
        {"layout_and_values", test_layout_and_values},
        {"determinism_and_order", test_determinism_and_order},
        {"multi_layer_bands", test_multi_layer_bands},
        {"precision_contract", test_precision_contract},
        {"parameters_key_rejected", test_parameters_key_rejected},
        {"rejections", test_rejections},
#ifdef FBS_ADAPTER_STUB_REGION
        {"stub_seams", test_stub_seams},
#endif
    };
    for (const auto& t : tests) {
        const int before = failures;
        try {
            t.second();
        } catch (const std::exception& e) {
            ++failures;
            std::cerr << t.first << ": unexpected exception: " << e.what() << "\n";
        }
        std::cout << (failures == before ? "PASS " : "FAIL ") << t.first << "\n";
    }
#ifndef FBS_ADAPTER_STUB_REGION
    std::vector<Region> generated;
    {
        const int before = failures;
        test_generated(fixtures.empty() ? nullptr : &generated);
        std::cout << (failures == before ? "PASS " : "FAIL ") << "generated_regions\n";
    }
#endif

    if (!fixtures.empty()) {
        bool ok = true;
        try {
            ok &= write_file(fixtures + "/constructed.csv",
                             to_moocow_csv({constructed({-3, 2, -1}), neighbour({-2, 2, -1})}, map_id,
                                           "Adapter, \"constructed\" map"));
            ok &= write_file(fixtures + "/layers.csv",
                             to_moocow_csv({with_id_group(constructed({0, 0, 0}), "1111"),
                                            with_id_group(neighbour({0, 0, 1}), "1111")},
                                           "5f0c7d2e-9a41-4b8e-8d3c-2e7f61a0b9c5", "Two layers"));
#ifndef FBS_ADAPTER_STUB_REGION
            if (!generated.empty())
                ok &= write_file(fixtures + "/generated.csv",
                                 to_moocow_csv(generated, "5f0c7d2e-9a41-4b8e-8d3c-2e7f61a0b9c7", "Generated pair"));
#endif
        } catch (const std::exception& e) {
            std::cerr << "fixture export failed: " << e.what() << "\n";
            ok = false;
        }
        if (!ok) {
            ++failures;
            std::cerr << "could not write fixtures to " << fixtures << "\n";
        } else {
            std::cout << "fixtures written to " << fixtures << "\n";
        }
    }

    std::cout << checks << " checks, " << failures << " failures\n";
    return failures == 0 ? 0 : 1;
}
