// Existing-format adapter: fbs::caves::Region -> MOOCoW sectioned CSV.
//
// The output is the exact layout written by MOOCoW's CsvExport.ExportDataset and
// read by CsvImport.ImportDataset (MOOCoW checkout eb3318a). No section, column,
// schema field or record kind is added. Every rule that MOOCoW's CSV parser,
// domain constructors and CaveDatasetValidator would apply is checked here first,
// so a returned string imports and validates without any MOOCoW change, and no
// supported value is silently altered or dropped by the importer. Anything the
// existing format cannot carry losslessly is rejected with an explanation.
//
// Coordinates: region-local meters (X/Y horizontal, +Z up) become MOOCoW world
// meters by adding Region::origin. Stitched tunnels are local to the first region
// passed to stitch() and are offset by that region's origin. See
// FinalBuildCaves/docs/MOOCOW-ADAPTER.md for the complete field mapping.

#include "fbs/caves.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace fbs::caves {
namespace {

// MOOCoW CaveDatasetValidator.GeometryTolerance and UnitVector3 tolerance.
constexpr double geometry_tolerance = 1e-9;
constexpr double unit_tolerance = 1e-12;
constexpr std::size_t max_reported_issues = 25;

[[noreturn]] void fail(ErrorCode code, const std::string& message) {
    throw Error(code, "to_moocow_csv: " + message);
}

// Products are forced through memory so the compiler cannot contract a*b+c
// into a fused multiply-add. MOOCoW's .NET arithmetic is unfused; validation
// decisions made here must be bit-identical to the ones the importer makes.
double mul(double a, double b) {
    volatile double r = a * b;
    return r;
}

// ---------------------------------------------------------------- text rules

// Decodes strict UTF-8 (no overlongs, surrogates or values above U+10FFFF).
bool decode_utf8(const std::string& s, std::vector<char32_t>& out) {
    out.clear();
    std::size_t i = 0;
    while (i < s.size()) {
        const auto c = static_cast<unsigned char>(s[i]);
        char32_t cp;
        std::size_t n;
        if (c < 0x80) { cp = c; n = 1; }
        else if (c >= 0xC2 && c <= 0xDF) { cp = c & 0x1F; n = 2; }
        else if (c >= 0xE0 && c <= 0xEF) { cp = c & 0x0F; n = 3; }
        else if (c >= 0xF0 && c <= 0xF4) { cp = c & 0x07; n = 4; }
        else return false;
        if (i + n > s.size()) return false;
        for (std::size_t k = 1; k < n; ++k) {
            const auto t = static_cast<unsigned char>(s[i + k]);
            if ((t & 0xC0) != 0x80) return false;
            cp = (cp << 6) | (t & 0x3F);
        }
        if ((n == 3 && cp < 0x800) || (n == 4 && (cp < 0x10000 || cp > 0x10FFFF)) ||
            (cp >= 0xD800 && cp <= 0xDFFF))
            return false;
        out.push_back(cp);
        i += n;
    }
    return true;
}

// .NET char.IsWhiteSpace for every code point (used by string.Trim and
// ArgumentException.ThrowIfNullOrWhiteSpace in DomainGuard.RequiredText).
bool dotnet_whitespace(char32_t c) {
    return (c >= 0x09 && c <= 0x0D) || c == 0x20 || c == 0x85 || c == 0xA0 || c == 0x1680 ||
           (c >= 0x2000 && c <= 0x200A) || c == 0x2028 || c == 0x2029 || c == 0x202F ||
           c == 0x205F || c == 0x3000;
}

void append_utf16(char32_t c, std::u16string& out) {
    if (c < 0x10000) {
        out.push_back(static_cast<char16_t>(c));
    } else {
        c -= 0x10000;
        out.push_back(static_cast<char16_t>(0xD800 + (c >> 10)));
        out.push_back(static_cast<char16_t>(0xDC00 + (c & 0x3FF)));
    }
}

std::u16string to_utf16(const std::string& s) {
    std::vector<char32_t> cps;
    decode_utf8(s, cps);  // callers validated already
    std::u16string out;
    for (auto c : cps) append_utf16(c, out);
    return out;
}

// Text stored in a plain CSV cell. CsvImport splits the file into lines before
// honouring quotes, so CR/LF can never survive; DomainGuard.RequiredText trims
// and rejects blank values, so surrounding whitespace would be silently lost.
// Tags and generation keys travel inside JSON, where line breaks are escaped,
// so only plain cells need the line-break rule.
void check_required_text(const std::string& value, const std::string& what, bool plain_cell = true) {
    std::vector<char32_t> cps;
    if (!decode_utf8(value, cps)) fail(ErrorCode::invalid_input, what + " is not valid UTF-8.");
    if (cps.empty()) fail(ErrorCode::invalid_input, what + " must not be empty.");
    if (dotnet_whitespace(cps.front()) || dotnet_whitespace(cps.back()))
        fail(ErrorCode::constraint,
             what + " has leading or trailing whitespace; MOOCoW trims it on import, so it cannot be preserved.");
    for (auto c : cps)
        if (plain_cell && (c == U'\r' || c == U'\n'))
            fail(ErrorCode::constraint,
                 what + " contains a line break; MOOCoW's CSV importer splits lines before parsing quotes.");
}

// Only the ASCII range is folded. .NET OrdinalIgnoreCase also folds non-ASCII
// letters; see tag rules (non-ASCII tags are rejected) and the documented gap
// for codes (MOOCoW's validator still rejects such collisions loudly).
std::string ascii_upper(const std::string& s) {
    std::string r = s;
    for (auto& ch : r)
        if (ch >= 'a' && ch <= 'z') ch = static_cast<char>(ch - 'a' + 'A');
    return r;
}

// ---------------------------------------------------------------- formatting

// CsvText.Number: double.ToString("G17", InvariantCulture). Same significant
// digits and fixed/scientific switch as C's %.17g; .NET writes an upper-case E.
std::string number(double v) {
    if (!std::isfinite(v)) fail(ErrorCode::invalid_input, "non-finite number reached the writer.");
    char buffer[64];
    const int n = std::snprintf(buffer, sizeof buffer, "%.17g", v);
    std::string s(buffer, static_cast<std::size_t>(n));
    for (auto& ch : s)
        if (ch == 'e') ch = 'E';
    return s;
}

// CsvText.Escape.
std::string csv(const std::string& value) {
    if (value.empty()) return {};
    if (value.find_first_of(",\"\n\r") == std::string::npos) return value;
    std::string out = "\"";
    for (char ch : value) {
        if (ch == '"') out += "\"\"";
        else out += ch;
    }
    out += '"';
    return out;
}

// System.Text.Json string literal with the default JavaScriptEncoder, as used
// by CsvText.FormatTags / FormatGenerationParameters. Only printable ASCII that
// is not HTML-sensitive stays literal; everything else is \uXXXX (upper-case
// hex, UTF-16 surrogate pairs) except the short escapes \b \t \n \f \r \\.
void append_hex4(std::string& out, unsigned v) {
    static const char digits[] = "0123456789ABCDEF";
    out += "\\u";
    for (int shift = 12; shift >= 0; shift -= 4) out += digits[(v >> shift) & 0xF];
}

std::string json_string(const std::string& value) {
    std::string out = "\"";
    for (char16_t u : to_utf16(value)) {
        switch (u) {
        case u'\b': out += "\\b"; continue;
        case u'\t': out += "\\t"; continue;
        case u'\n': out += "\\n"; continue;
        case u'\f': out += "\\f"; continue;
        case u'\r': out += "\\r"; continue;
        case u'\\': out += "\\\\"; continue;
        default: break;
        }
        const bool literal = u >= 0x20 && u <= 0x7E && u != u'"' && u != u'&' && u != u'\'' &&
                             u != u'<' && u != u'>' && u != u'+' && u != u'`';
        if (literal) out += static_cast<char>(u);
        else append_hex4(out, u);
    }
    out += '"';
    return out;
}

bool ordinal_ignore_case_less(const std::string& a, const std::string& b) {
    return ascii_upper(a) < ascii_upper(b);
}

// Tags as MOOCoW stores them: trimmed, distinct ignoring case, sorted with
// OrdinalIgnoreCase. Order is not a stored attribute, so the canonical order is
// written. Anything MOOCoW would change or merge is rejected instead.
std::vector<std::string> canonical_tags(const std::vector<std::string>& tags, const std::string& owner) {
    std::vector<std::string> out;
    std::set<std::string> folded;
    for (const auto& tag : tags) {
        const std::string what = owner + " tag '" + tag + "'";
        check_required_text(tag, what, false);
        for (char ch : tag)
            if (static_cast<unsigned char>(ch) >= 0x80)
                fail(ErrorCode::constraint,
                     what + " contains non-ASCII text; MOOCoW deduplicates tags with .NET OrdinalIgnoreCase "
                            "and the adapter cannot prove it would not merge it with another tag.");
        if (!folded.insert(ascii_upper(tag)).second)
            fail(ErrorCode::constraint,
                 what + " duplicates another tag ignoring case; MOOCoW would silently merge them.");
        out.push_back(tag);
    }
    std::sort(out.begin(), out.end(), ordinal_ignore_case_less);
    return out;
}

std::string tags_json(const std::vector<std::string>& canonical) {
    std::string out = "[";
    for (std::size_t i = 0; i < canonical.size(); ++i) {
        if (i) out += ',';
        out += json_string(canonical[i]);
    }
    return out + "]";
}

// GenerationSettings keeps an ordinal (UTF-16) sorted copy with trimmed keys.
// A key named exactly "parameters" is refused: CsvText.ParseGenerationParameters
// and the Cavern/Tunnel repositories treat a root "parameters" property as a
// legacy wrapper, so MOOCoW either rejects the CSV or stores a record whose
// map can no longer be loaded (verified by adapter_acceptance.cs). No existing
// JSON shape carries that key through CSV import, SQLite storage and reload.
std::string generation_json(const Generation& generation, const std::string& owner) {
    std::vector<std::pair<std::u16string, const std::pair<const std::string, std::string>*>> sorted;
    for (const auto& entry : generation.parameters) {
        check_required_text(entry.first, owner + " generation parameter key '" + entry.first + "'", false);
        std::vector<char32_t> cps;
        if (!decode_utf8(entry.second, cps))
            fail(ErrorCode::invalid_input, owner + " generation parameter '" + entry.first + "' is not valid UTF-8.");
        if (entry.first == "parameters")
            fail(ErrorCode::constraint,
                 owner + " has a generation parameter named 'parameters'; MOOCoW reads a root 'parameters' "
                         "property as a legacy wrapper and cannot store and reload it. Rename the key.");
        sorted.emplace_back(to_utf16(entry.first), &entry);
    }
    std::sort(sorted.begin(), sorted.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });
    std::string body = "{";
    for (std::size_t i = 0; i < sorted.size(); ++i) {
        if (i) body += ',';
        body += json_string(sorted[i].second->first) + ':' + json_string(sorted[i].second->second);
    }
    return body + '}';
}

const char* authoring_text(Authoring a) {
    switch (a) {
    case Authoring::exact: return "exact";
    case Authoring::constrained_procedural: return "constrained_procedural";
    case Authoring::generated: return "generated";
    }
    fail(ErrorCode::invalid_input, "unknown authoring mode.");
}

const char* traversal_text(Traversal t) {
    switch (t) {
    case Traversal::walkable: return "walkable";
    case Traversal::difficult: return "difficult";
    case Traversal::hazardous: return "hazardous";
    case Traversal::restricted: return "restricted";
    case Traversal::impassable: return "impassable";
    }
    fail(ErrorCode::invalid_input, "unknown traversal class.");
}

const char* transition_text(Transition t) {
    switch (t) {
    case Transition::walkable: return "walkable";
    case Transition::slope: return "slope";
    case Transition::ramp: return "ramp";
    case Transition::stairs: return "stairs";
    case Transition::bridge: return "bridge";
    case Transition::climbable_cliff: return "climbable_cliff";
    case Transition::impassable_cliff: return "impassable_cliff";
    case Transition::one_way_drop: return "one_way_drop";
    case Transition::unconnected: return "unconnected";
    }
    fail(ErrorCode::invalid_input, "unknown transition kind.");
}

const char* port_state_text(PortState s) {
    switch (s) {
    case PortState::required: return "required";
    case PortState::optional: return "optional";
    case PortState::generated: return "generated";
    case PortState::sealed: return "sealed";
    }
    fail(ErrorCode::invalid_input, "unknown port state.");
}

// ---------------------------------------------------------------- identities

bool hex_digit(char c) {
    return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f') || (c >= 'A' && c <= 'F');
}

// Canonical 8-4-4-4-12 UUID, any hex case, not nil. Returned in the lower-case
// "D" form EntityId.ToString() writes; the 128-bit value is unchanged.
std::string uuid(const std::string& text, const std::string& what) {
    bool ok = text.size() == 36;
    bool nonzero = false;
    for (std::size_t i = 0; ok && i < text.size(); ++i) {
        if (i == 8 || i == 13 || i == 18 || i == 23) ok = text[i] == '-';
        else {
            ok = hex_digit(text[i]);
            nonzero = nonzero || text[i] != '0';
        }
    }
    if (!ok) fail(ErrorCode::invalid_input, what + " '" + text + "' is not a canonical UUID string.");
    if (!nonzero) fail(ErrorCode::invalid_input, what + " must not be the nil UUID.");
    std::string out = text;
    for (auto& ch : out)
        if (ch >= 'A' && ch <= 'F') ch = static_cast<char>(ch - 'A' + 'a');
    return out;
}

std::string optional_uuid(const std::string& text, const std::string& what) {
    return text.empty() ? std::string() : uuid(text, what);
}

// ---------------------------------------------------------------- numbers

void finite(double v, const std::string& what) {
    if (!std::isfinite(v)) fail(ErrorCode::invalid_input, what + " must be finite.");
}

void finite(const Vec3& v, const std::string& what) {
    finite(v.x, what + ".x");
    finite(v.y, what + ".y");
    finite(v.z, what + ".z");
}

double shifted(double local, double offset, const std::string& what) {
    finite(local, what);
    const double world = local + offset;
    if (!std::isfinite(world)) fail(ErrorCode::invalid_input, what + " overflows when offset to world coordinates.");
    return world;
}

Vec3 shifted(const Vec3& v, const Vec3& o, const std::string& what) {
    return {shifted(v.x, o.x, what + ".x"), shifted(v.y, o.y, what + ".y"), shifted(v.z, o.z, what + ".z")};
}

Vec2 shifted(const Vec2& v, const Vec3& o, const std::string& what) {
    return {shifted(v.x, o.x, what + ".x"), shifted(v.y, o.y, what + ".y")};
}

// ---------------------------------------------------------------- world model

struct WorldCavern { Cavern c; std::vector<std::string> tags; std::string gen; };
struct WorldFloor { Floor f; std::vector<std::string> tags; };
struct WorldTunnel { Tunnel t; std::vector<std::string> tags; std::string gen; };
struct Band { std::string id, code, name; double min_z, max_z; int order; };

struct World {
    std::string map_id, map_code, map_name;
    Vec3 min{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(),
             std::numeric_limits<double>::infinity()};
    Vec3 max{-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(),
             -std::numeric_limits<double>::infinity()};
    std::vector<Band> bands;
    std::vector<WorldCavern> caverns;
    std::vector<WorldFloor> floors;
    std::vector<Cliff> cliffs;
    std::vector<Port> ports;
    std::vector<WorldTunnel> tunnels;

    void include(double x, double y, double z) {
        min.x = std::min(min.x, x); min.y = std::min(min.y, y); min.z = std::min(min.z, z);
        max.x = std::max(max.x, x); max.y = std::max(max.y, y); max.z = std::max(max.z, z);
    }
    void include_xy(double x, double y) {
        min.x = std::min(min.x, x); min.y = std::min(min.y, y);
        max.x = std::max(max.x, x); max.y = std::max(max.y, y);
    }
    void include_z(double z) {
        min.z = std::min(min.z, z);
        max.z = std::max(max.z, z);
    }
};

std::string label(const char* kind, const std::string& id, const std::string& code) {
    return std::string(kind) + " " + id + (code.empty() ? "" : " ('" + code + "')");
}

Cavern world_cavern(const Cavern& in, const Vec3& o) {
    const std::string what = label("cavern", in.id, in.code);
    Cavern c = in;
    c.id = uuid(in.id, what + " id");
    check_required_text(in.code, what + " code");
    check_required_text(in.name, what + " name");
    c.center = shifted(in.center, o, what + " center");
    c.bounds.minimum = shifted(in.bounds.minimum, o, what + " bounds.minimum");
    c.bounds.maximum = shifted(in.bounds.maximum, o, what + " bounds.maximum");
    finite(in.reserved_clearance, what + " reserved_clearance");
    for (std::size_t i = 0; i < in.boundary.size(); ++i)
        c.boundary[i] = shifted(in.boundary[i], o, what + " outline point");
    return c;
}

Floor world_floor(const Floor& in, const Vec3& o) {
    const std::string what = label("floor", in.id, in.code);
    Floor f = in;
    f.id = uuid(in.id, what + " id");
    f.cavern_id = uuid(in.cavern_id, what + " cavern_id");
    check_required_text(in.code, what + " code");
    check_required_text(in.name, what + " name");
    check_required_text(in.material, what + " material");
    for (std::size_t i = 0; i < in.boundary.size(); ++i)
        f.boundary[i] = shifted(in.boundary[i], o, what + " boundary point");
    // Slopes are gradients anchored at the floor's own first boundary vertex,
    // so a translation changes only base_z (see MOOCOW-ADAPTER.md).
    f.base_z = shifted(in.base_z, o.z, what + " base_z");
    finite(in.slope_x, what + " slope_x");
    finite(in.slope_y, what + " slope_y");
    finite(in.variation_amplitude, what + " variation_amplitude");
    return f;
}

Cliff world_cliff(const Cliff& in, const Vec3& o) {
    const std::string what = label("cliff", in.id, in.code);
    Cliff c = in;
    c.id = uuid(in.id, what + " id");
    c.cavern_id = uuid(in.cavern_id, what + " cavern_id");
    c.from_floor_id = uuid(in.from_floor_id, what + " from_floor_id");
    c.to_floor_id = optional_uuid(in.to_floor_id, what + " to_floor_id");
    check_required_text(in.code, what + " code");
    for (std::size_t i = 0; i < in.boundary.size(); ++i)
        c.boundary[i] = shifted(in.boundary[i], o, what + " point");
    finite(in.height, what + " height");
    return c;
}

Port world_port(const Port& in, const Vec3& o) {
    const std::string what = label("port", in.id, in.code);
    Port p = in;
    p.id = uuid(in.id, what + " id");
    p.cavern_id = uuid(in.cavern_id, what + " cavern_id");
    p.floor_id = optional_uuid(in.floor_id, what + " floor_id");
    check_required_text(in.code, what + " code");
    check_required_text(in.name, what + " name");
    p.position = shifted(in.position, o, what + " position");
    finite(in.facing, what + " facing");
    finite(in.width, what + " width");
    finite(in.height, what + " height");
    return p;
}

Tunnel world_tunnel(const Tunnel& in, const Vec3& o) {
    const std::string what = label("tunnel", in.id, in.code);
    Tunnel t = in;
    t.id = uuid(in.id, what + " id");
    t.start_port_id = uuid(in.start_port_id, what + " start_port_id");
    t.end_port_id = uuid(in.end_port_id, what + " end_port_id");
    check_required_text(in.code, what + " code");
    check_required_text(in.name, what + " name");
    for (std::size_t i = 0; i < in.centerline.size(); ++i) {
        t.centerline[i].center = shifted(in.centerline[i].center, o, what + " profile center");
        finite(in.centerline[i].width, what + " profile width");
        finite(in.centerline[i].height, what + " profile height");
    }
    finite(in.maximum_slope_degrees, what + " maximum_slope_degrees");
    finite(in.minimum_clearance, what + " minimum_clearance");
    return t;
}

// ---------------------------------------------------------------- regions

bool same_settings(const Region& a, const Region& b) {
    return a.version == b.version && a.algorithm == b.algorithm && a.world_seed == b.world_seed &&
           a.settings_id == b.settings_id && a.size.x == b.size.x && a.size.y == b.size.y &&
           a.size.z == b.size.z;
}

bool key_less(const RegionKey& a, const RegionKey& b) {
    if (a.z != b.z) return a.z < b.z;
    if (a.y != b.y) return a.y < b.y;
    return a.x < b.x;
}

std::string key_text(const RegionKey& k) {
    return "(" + std::to_string(k.x) + "," + std::to_string(k.y) + "," + std::to_string(k.z) + ")";
}

// Origins must lie on one common lattice: origin - key*size is the same for
// every region, so shared faces coincide in world space.
void check_lattice(const Region& base, const Region& r) {
    const double dk[3] = {static_cast<double>(r.key.x) - static_cast<double>(base.key.x),
                          static_cast<double>(r.key.y) - static_cast<double>(base.key.y),
                          static_cast<double>(r.key.z) - static_cast<double>(base.key.z)};
    const double size[3] = {base.size.x, base.size.y, base.size.z};
    const double o0[3] = {base.origin.x, base.origin.y, base.origin.z};
    const double o1[3] = {r.origin.x, r.origin.y, r.origin.z};
    for (int a = 0; a < 3; ++a) {
        const double expected = dk[a] * size[a];
        const double actual = o1[a] - o0[a];
        const double scale = std::max({1.0, std::fabs(o0[a]), std::fabs(o1[a]), std::fabs(expected)});
        if (!(std::fabs(actual - expected) <= 1e-9 + 8 * std::numeric_limits<double>::epsilon() * scale))
            fail(ErrorCode::incompatible_boundary,
                 "region " + key_text(r.key) + " origin is not on the lattice of region " + key_text(base.key) +
                     " (origin - key*size differs).");
    }
}

// ---------------------------------------------------------------- validation

struct P2 { double x, y; };

double orientation(P2 a, P2 b, P2 c) {
    return mul(b.x - a.x, c.y - a.y) - mul(b.y - a.y, c.x - a.x);
}

bool boxes_overlap(P2 a, P2 b, P2 c, P2 d) {
    return std::max(std::min(a.x, b.x), std::min(c.x, d.x)) <=
               std::min(std::max(a.x, b.x), std::max(c.x, d.x)) + geometry_tolerance &&
           std::max(std::min(a.y, b.y), std::min(c.y, d.y)) <=
               std::min(std::max(a.y, b.y), std::max(c.y, d.y)) + geometry_tolerance;
}

bool segments_intersect(P2 a, P2 b, P2 c, P2 d) {
    return mul(orientation(a, b, c), orientation(a, b, d)) <= geometry_tolerance &&
           mul(orientation(c, d, a), orientation(c, d, b)) <= geometry_tolerance && boxes_overlap(a, b, c, d);
}

bool on_segment(P2 a, P2 b, P2 p) {
    return std::fabs(orientation(a, b, p)) <= geometry_tolerance && boxes_overlap(a, b, p, p);
}

bool polygon_contains(const Polygon& polygon, P2 p) {
    bool inside = false;
    const std::size_t n = polygon.size();
    for (std::size_t i = 0, j = n - 1; i < n; j = i++) {
        const P2 a{polygon[i].x, polygon[i].y};
        const P2 b{polygon[j].x, polygon[j].y};
        if (on_segment(a, b, p)) return true;
        if ((a.y > p.y) != (b.y > p.y) && p.x < (b.x - a.x) * (p.y - a.y) / (b.y - a.y) + a.x) inside = !inside;
    }
    return inside;
}

double distance(const Vec3& a, const Vec3& b) {
    return std::sqrt(std::pow(a.x - b.x, 2) + std::pow(a.y - b.y, 2) + std::pow(a.z - b.z, 2));
}

bool inside_xy(const Bounds& b, const Vec2& p) {
    return p.x >= b.minimum.x && p.x <= b.maximum.x && p.y >= b.minimum.y && p.y <= b.maximum.y;
}

bool inside_3d(const Bounds& b, const Vec3& p) {
    return p.x >= b.minimum.x && p.x <= b.maximum.x && p.y >= b.minimum.y && p.y <= b.maximum.y &&
           p.z >= b.minimum.z && p.z <= b.maximum.z;
}

bool has_isolated_tag(const std::vector<std::string>& tags) {
    for (const auto& t : tags)
        if (ascii_upper(t) == "ISOLATED") return true;
    return false;
}

struct Issues {
    std::vector<std::pair<ErrorCode, std::string>> list;
    void add(ErrorCode code, const std::string& message) { list.emplace_back(code, message); }
    void raise() const {
        if (list.empty()) return;
        std::string text = std::to_string(list.size()) + " MOOCoW rule violation(s):";
        for (std::size_t i = 0; i < list.size() && i < max_reported_issues; ++i) text += "\n  - " + list[i].second;
        if (list.size() > max_reported_issues) text += "\n  - ...";
        fail(list.front().first, text);
    }
};

// Distinct-vertex and ordering rules of PlanarPolygon plus the validator's
// consecutive-duplicate and self-intersection checks.
void check_polygon(const Polygon& points, const std::string& what, Issues& issues) {
    std::set<std::pair<double, double>> distinct;
    for (const auto& p : points) distinct.emplace(p.x == 0 ? 0.0 : p.x, p.y == 0 ? 0.0 : p.y);
    if (points.size() < 3 || distinct.size() < 3) {
        issues.add(ErrorCode::invalid_input, what + " needs at least three distinct vertices.");
        return;
    }
    for (std::size_t i = 0; i + 1 < points.size(); ++i)
        if (points[i].x == points[i + 1].x && points[i].y == points[i + 1].y) {
            issues.add(ErrorCode::invalid_input, what + " contains consecutive duplicate points.");
            break;
        }
    const std::size_t n = points.size();
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = i + 1; j < n; ++j) {
            if (j == i + 1 || (i == 0 && j == n - 1)) continue;
            if (segments_intersect({points[i].x, points[i].y}, {points[(i + 1) % n].x, points[(i + 1) % n].y},
                                   {points[j].x, points[j].y}, {points[(j + 1) % n].x, points[(j + 1) % n].y})) {
                issues.add(ErrorCode::invalid_input, what + " edges intersect.");
                return;
            }
        }
}

void check_identity_and_codes(const World& w, Issues& issues) {
    std::unordered_map<std::string, int> ids;
    // Folded "scope\0code" -> (count, first original scope/code for the message).
    std::unordered_map<std::string, std::pair<int, std::pair<std::string, std::string>>> codes;
    auto id = [&](const std::string& value) { ++ids[value]; };
    auto code = [&](const std::string& scope, const std::string& value) {
        auto& entry = codes[ascii_upper(scope + '\0' + value)];
        if (entry.first++ == 0) entry.second = {scope, value};
    };
    id(w.map_id);
    for (const auto& b : w.bands) { id(b.id); code("map:" + w.map_id, b.code); }
    for (const auto& c : w.caverns) { id(c.c.id); code("map:" + w.map_id, c.c.code); }
    for (const auto& f : w.floors) { id(f.f.id); code("cavern:" + f.f.cavern_id, f.f.code); }
    for (const auto& c : w.cliffs) { id(c.id); code("cavern:" + c.cavern_id, c.code); }
    for (const auto& p : w.ports) { id(p.id); code("cavern:" + p.cavern_id, p.code); }
    for (const auto& t : w.tunnels) { id(t.t.id); code("map:" + w.map_id, t.t.code); }
    std::vector<std::string> dup_ids, dup_codes;
    for (const auto& e : ids)
        if (e.second > 1) dup_ids.push_back(e.first + " is used " + std::to_string(e.second) + " times");
    for (const auto& e : codes)
        if (e.second.first > 1)
            dup_codes.push_back("code '" + e.second.second.second + "' (case-insensitive) is duplicated in " +
                                e.second.second.first);
    std::sort(dup_ids.begin(), dup_ids.end());
    std::sort(dup_codes.begin(), dup_codes.end());
    for (const auto& m : dup_ids) issues.add(ErrorCode::invalid_input, "duplicate entity ID: " + m + ".");
    for (const auto& m : dup_codes) issues.add(ErrorCode::constraint, m + ".");
}

void check_records(const World& w, Issues& issues) {
    std::unordered_map<std::string, const WorldCavern*> caverns;
    std::unordered_map<std::string, const WorldFloor*> floors;
    std::unordered_map<std::string, const Port*> ports;
    for (const auto& c : w.caverns) caverns.emplace(c.c.id, &c);
    for (const auto& f : w.floors) floors.emplace(f.f.id, &f);
    for (const auto& p : w.ports) ports.emplace(p.id, &p);

    for (const auto& wc : w.caverns) {
        const Cavern& c = wc.c;
        const std::string what = label("cavern", c.id, c.code);
        const Bounds& b = c.bounds;
        if (b.minimum.x > b.maximum.x || b.minimum.y > b.maximum.y || b.minimum.z > b.maximum.z)
            issues.add(ErrorCode::invalid_input, what + " bounds minimum exceeds maximum.");
        else if (!std::isfinite(b.maximum.x - b.minimum.x) || !std::isfinite(b.maximum.y - b.minimum.y) ||
                 !std::isfinite(b.maximum.z - b.minimum.z))
            issues.add(ErrorCode::invalid_input, what + " bounds extent is not finite.");
        else if (!inside_3d(b, c.center))
            issues.add(ErrorCode::invalid_input, what + " center is outside its bounds.");
        if (c.reserved_clearance < 0) issues.add(ErrorCode::invalid_input, what + " reserved_clearance is negative.");
        if (!c.boundary.empty()) {
            check_polygon(c.boundary, what + " outline", issues);
            for (const auto& p : c.boundary)
                if (!inside_xy(b, p)) {
                    issues.add(ErrorCode::invalid_input, what + " outline is outside its bounds.");
                    break;
                }
        } else if (c.authoring == Authoring::exact) {
            issues.add(ErrorCode::invalid_input, what + " is exact but has no outline.");
        }
    }

    for (const auto& wf : w.floors) {
        const Floor& f = wf.f;
        const std::string what = label("floor", f.id, f.code);
        if (f.variation_amplitude < 0) issues.add(ErrorCode::invalid_input, what + " variation_amplitude is negative.");
        check_polygon(f.boundary, what + " boundary", issues);
        const auto owner = caverns.find(f.cavern_id);
        if (owner == caverns.end()) {
            issues.add(ErrorCode::invalid_input, what + " references a missing cavern.");
            continue;
        }
        for (const auto& p : f.boundary)
            if (!inside_xy(owner->second->c.bounds, p)) {
                issues.add(ErrorCode::invalid_input, what + " is outside its cavern bounds.");
                break;
            }
    }

    for (const auto& c : w.cliffs) {
        const std::string what = label("cliff", c.id, c.code);
        if (caverns.find(c.cavern_id) == caverns.end())
            issues.add(ErrorCode::invalid_input, what + " references a missing cavern.");
        const auto from = floors.find(c.from_floor_id);
        if (from == floors.end()) issues.add(ErrorCode::invalid_input, what + " references a missing source floor.");
        else if (from->second->f.cavern_id != c.cavern_id)
            issues.add(ErrorCode::invalid_input, what + " source floor belongs to another cavern.");
        if (!c.to_floor_id.empty()) {
            const auto to = floors.find(c.to_floor_id);
            if (to == floors.end()) issues.add(ErrorCode::invalid_input, what + " references a missing destination floor.");
            else if (to->second->f.cavern_id != c.cavern_id)
                issues.add(ErrorCode::invalid_input, what + " destination floor belongs to another cavern.");
        }
        if (c.to_floor_id == c.from_floor_id)
            issues.add(ErrorCode::invalid_input, what + " connects a floor to itself.");
        if (c.boundary.size() < 2) issues.add(ErrorCode::invalid_input, what + " needs at least two points.");
        for (std::size_t i = 0; i + 1 < c.boundary.size(); ++i)
            if (c.boundary[i].x == c.boundary[i + 1].x && c.boundary[i].y == c.boundary[i + 1].y) {
                issues.add(ErrorCode::invalid_input, what + " contains consecutive duplicate points.");
                break;
            }
        if (c.height < 0) issues.add(ErrorCode::invalid_input, what + " height is negative.");
        if (c.transition == Transition::one_way_drop && c.to_floor_id.empty())
            issues.add(ErrorCode::constraint, what + " is a one-way drop without a destination floor.");
        const bool needs_height = c.transition == Transition::slope || c.transition == Transition::ramp ||
                                  c.transition == Transition::stairs || c.transition == Transition::climbable_cliff ||
                                  c.transition == Transition::impassable_cliff ||
                                  c.transition == Transition::one_way_drop;
        if (needs_height && c.height <= geometry_tolerance)
            issues.add(ErrorCode::constraint, what + " transition requires a positive height.");
    }

    std::unordered_set<std::string> used_ports;
    for (const auto& t : w.tunnels) {
        used_ports.insert(t.t.start_port_id);
        used_ports.insert(t.t.end_port_id);
    }

    for (const auto& p : w.ports) {
        const std::string what = label("port", p.id, p.code);
        if (!(p.width > 0) || !(p.height > 0)) issues.add(ErrorCode::invalid_input, what + " width and height must be positive.");
        const double m = std::sqrt(mul(p.facing.x, p.facing.x) + mul(p.facing.y, p.facing.y) + mul(p.facing.z, p.facing.z));
        if (!(m > unit_tolerance)) {
            issues.add(ErrorCode::invalid_input, what + " facing has zero magnitude.");
        } else if (p.facing.x / m != p.facing.x || p.facing.y / m != p.facing.y || p.facing.z / m != p.facing.z) {
            issues.add(ErrorCode::constraint,
                       what + " facing is not an exactly normalized unit vector; MOOCoW's UnitVector3 would "
                              "rescale it on import.");
        }
        const auto owner = caverns.find(p.cavern_id);
        if (owner == caverns.end()) {
            issues.add(ErrorCode::invalid_input, what + " references a missing cavern.");
        } else if (!inside_3d(owner->second->c.bounds, p.position)) {
            issues.add(ErrorCode::constraint, what + " is outside its cavern bounds.");
        }
        if (!p.floor_id.empty()) {
            const auto floor = floors.find(p.floor_id);
            if (floor == floors.end()) issues.add(ErrorCode::invalid_input, what + " references a missing floor.");
            else {
                if (floor->second->f.cavern_id != p.cavern_id)
                    issues.add(ErrorCode::invalid_input, what + " floor belongs to another cavern.");
                if (floor->second->f.boundary.size() >= 3 &&
                    !polygon_contains(floor->second->f.boundary, {p.position.x, p.position.y}))
                    issues.add(ErrorCode::constraint, what + " is outside its floor region.");
            }
        }
        const bool used = used_ports.count(p.id) != 0;
        if (p.state == PortState::required && !used)
            issues.add(ErrorCode::constraint, what + " is required but no tunnel uses it.");
        if (p.state == PortState::sealed && used)
            issues.add(ErrorCode::constraint, what + " is sealed but a tunnel uses it.");
    }

    for (const auto& wt : w.tunnels) {
        const Tunnel& t = wt.t;
        const std::string what = label("tunnel", t.id, t.code);
        if (t.start_port_id == t.end_port_id) issues.add(ErrorCode::invalid_input, what + " must connect two different ports.");
        if (!(t.maximum_slope_degrees >= 0 && t.maximum_slope_degrees < 90))
            issues.add(ErrorCode::invalid_input, what + " maximum slope must be in [0, 90) degrees.");
        if (!(t.minimum_clearance > 0)) issues.add(ErrorCode::invalid_input, what + " minimum clearance must be positive.");
        if (t.centerline.size() < 2) {
            issues.add(ErrorCode::invalid_input, what + " needs at least two centerline profiles.");
            continue;
        }
        for (const auto& profile : t.centerline) {
            if (!(profile.width > 0) || !(profile.height > 0)) {
                issues.add(ErrorCode::invalid_input, what + " profile width and height must be positive.");
                break;
            }
            if (profile.height < t.minimum_clearance) {
                issues.add(ErrorCode::constraint, what + " profile height is below its minimum clearance.");
                break;
            }
        }
        const auto start = ports.find(t.start_port_id);
        const auto end = ports.find(t.end_port_id);
        if (start == ports.end()) issues.add(ErrorCode::invalid_input, what + " start references a missing port.");
        if (end == ports.end()) issues.add(ErrorCode::invalid_input, what + " end references a missing port.");
        if (start != ports.end() && end != ports.end()) {
            const Port& s = *start->second;
            const Port& e = *end->second;
            if (distance(t.centerline.front().center, s.position) > std::max(s.width, s.height))
                issues.add(ErrorCode::constraint, what + " misses its start port.");
            if (distance(t.centerline.back().center, e.position) > std::max(e.width, e.height))
                issues.add(ErrorCode::constraint, what + " misses its end port.");
            for (std::size_t i = 1; i < t.centerline.size(); ++i) {
                const Vec3& a = t.centerline[i - 1].center;
                const Vec3& b = t.centerline[i].center;
                const double horizontal = std::sqrt(std::pow(b.x - a.x, 2) + std::pow(b.y - a.y, 2));
                const double slope = horizontal <= geometry_tolerance
                                         ? 90
                                         : std::atan2(std::fabs(b.z - a.z), horizontal) * 180 / 3.14159265358979323846;
                if (slope > t.maximum_slope_degrees + geometry_tolerance) {
                    issues.add(ErrorCode::constraint, what + " segment " + std::to_string(i) +
                                                          " exceeds its maximum slope in world coordinates.");
                    break;
                }
            }
        }
    }

    // Connectivity: every cavern needs a tunnel-used port or an "isolated" tag,
    // and floors inside one cavern must form one component through connecting
    // cliffs unless the members of the smaller components are tagged isolated.
    std::unordered_set<std::string> linked;
    for (const auto& id : used_ports) {
        const auto p = ports.find(id);
        if (p != ports.end()) linked.insert(p->second->cavern_id);
    }
    for (const auto& c : w.caverns)
        if (!linked.count(c.c.id) && !has_isolated_tag(c.tags))
            issues.add(ErrorCode::constraint,
                       label("cavern", c.c.id, c.c.code) + " has no connected ports and is not tagged isolated.");

    std::map<std::string, std::vector<const WorldFloor*>> by_cavern;
    for (const auto& f : w.floors) by_cavern[f.f.cavern_id].push_back(&f);
    for (const auto& group : by_cavern) {
        if (group.second.size() < 2) continue;
        std::unordered_map<std::string, std::vector<std::string>> edges;
        for (const auto* f : group.second) edges[f->f.id];
        bool broken = false;
        for (const auto& c : w.cliffs) {
            if (c.cavern_id != group.first || c.to_floor_id.empty() || c.transition == Transition::unconnected ||
                c.transition == Transition::impassable_cliff)
                continue;
            if (!edges.count(c.from_floor_id) || !edges.count(c.to_floor_id)) { broken = true; continue; }
            edges[c.from_floor_id].push_back(c.to_floor_id);
            edges[c.to_floor_id].push_back(c.from_floor_id);
        }
        if (broken) continue;  // already reported as a reference error
        std::unordered_set<std::string> remaining;
        for (const auto* f : group.second) remaining.insert(f->f.id);
        for (const auto* seed : group.second) {
            if (!remaining.count(seed->f.id)) continue;
            std::vector<std::string> component{seed->f.id}, stack{seed->f.id};
            remaining.erase(seed->f.id);
            while (!stack.empty()) {
                const std::string id = stack.back();
                stack.pop_back();
                for (const auto& next : edges[id])
                    if (remaining.erase(next)) {
                        component.push_back(next);
                        stack.push_back(next);
                    }
            }
            if (component.size() >= group.second.size()) continue;
            for (const auto& id : component) {
                const auto* f = floors.at(id);
                if (!has_isolated_tag(f->tags))
                    issues.add(ErrorCode::constraint, label("floor", f->f.id, f->f.code) +
                                                          " is in an isolated floor component without an isolated tag.");
            }
        }
    }
}

// ---------------------------------------------------------------- emit

std::string emit(const World& w) {
    std::string s;
    s.reserve(4096);
    auto line = [&](const std::string& text) { s += text; s += '\n'; };

    line("=== MAP ===");
    line("id,code,name,min_x,min_y,min_z,max_x,max_y,max_z");
    line(w.map_id + ',' + csv(w.map_code) + ',' + csv(w.map_name) + ',' + number(w.min.x) + ',' + number(w.min.y) +
         ',' + number(w.min.z) + ',' + number(w.max.x) + ',' + number(w.max.y) + ',' + number(w.max.z));

    line("");
    line("=== ELEVATION_BANDS ===");
    line("id,map_id,code,name,min_z,max_z,display_order");
    for (const auto& b : w.bands)
        line(b.id + ',' + w.map_id + ',' + csv(b.code) + ',' + csv(b.name) + ',' + number(b.min_z) + ',' +
             number(b.max_z) + ',' + std::to_string(b.order));

    line("");
    line("=== CAVERNS ===");
    line("id,map_id,code,name,center_x,center_y,center_z,min_x,min_y,min_z,max_x,max_y,max_z,"
         "authoring_mode,is_locked,reserved_clearance_m,generation_seed,tags,generation_json");
    for (const auto& wc : w.caverns) {
        const Cavern& c = wc.c;
        line(c.id + ',' + w.map_id + ',' + csv(c.code) + ',' + csv(c.name) + ',' + number(c.center.x) + ',' +
             number(c.center.y) + ',' + number(c.center.z) + ',' + number(c.bounds.minimum.x) + ',' +
             number(c.bounds.minimum.y) + ',' + number(c.bounds.minimum.z) + ',' + number(c.bounds.maximum.x) + ',' +
             number(c.bounds.maximum.y) + ',' + number(c.bounds.maximum.z) + ',' + authoring_text(c.authoring) + ',' +
             (c.locked ? "1" : "0") + ',' + number(c.reserved_clearance) + ',' + std::to_string(c.generation.seed) +
             ',' + csv(tags_json(wc.tags)) + ',' + csv(wc.gen));
    }

    line("");
    line("=== CAVERN_OUTLINE_POINTS ===");
    line("cavern_id,point_index,x,y");
    for (const auto& wc : w.caverns)
        for (std::size_t i = 0; i < wc.c.boundary.size(); ++i)
            line(wc.c.id + ',' + std::to_string(i) + ',' + number(wc.c.boundary[i].x) + ',' + number(wc.c.boundary[i].y));

    line("");
    line("=== FLOOR_REGIONS ===");
    line("id,cavern_id,code,name,base_z,slope_x,slope_y,variation_amplitude_m,variation_seed,"
         "material_code,traversal_class,tags");
    for (const auto& wf : w.floors) {
        const Floor& f = wf.f;
        line(f.id + ',' + f.cavern_id + ',' + csv(f.code) + ',' + csv(f.name) + ',' + number(f.base_z) + ',' +
             number(f.slope_x) + ',' + number(f.slope_y) + ',' + number(f.variation_amplitude) + ',' +
             std::to_string(f.variation_seed) + ',' + csv(f.material) + ',' + traversal_text(f.traversal) + ',' +
             csv(tags_json(wf.tags)));
    }

    line("");
    line("=== FLOOR_REGION_POINTS ===");
    line("floor_region_id,point_index,x,y");
    for (const auto& wf : w.floors)
        for (std::size_t i = 0; i < wf.f.boundary.size(); ++i)
            line(wf.f.id + ',' + std::to_string(i) + ',' + number(wf.f.boundary[i].x) + ',' + number(wf.f.boundary[i].y));

    line("");
    line("=== CLIFF_EDGES ===");
    line("id,cavern_id,from_floor_region_id,to_floor_region_id,code,transition_kind,height_m");
    for (const auto& c : w.cliffs)
        line(c.id + ',' + c.cavern_id + ',' + c.from_floor_id + ',' + c.to_floor_id + ',' + csv(c.code) + ',' +
             transition_text(c.transition) + ',' + number(c.height));

    line("");
    line("=== CLIFF_EDGE_POINTS ===");
    line("cliff_edge_id,point_index,x,y");
    for (const auto& c : w.cliffs)
        for (std::size_t i = 0; i < c.boundary.size(); ++i)
            line(c.id + ',' + std::to_string(i) + ',' + number(c.boundary[i].x) + ',' + number(c.boundary[i].y));

    line("");
    line("=== CAVERN_PORTS ===");
    line("id,cavern_id,floor_region_id,code,name,position_x,position_y,position_z,"
         "facing_x,facing_y,facing_z,width_m,height_m,port_state,traversal_class");
    for (const auto& p : w.ports)
        line(p.id + ',' + p.cavern_id + ',' + p.floor_id + ',' + csv(p.code) + ',' + csv(p.name) + ',' +
             number(p.position.x) + ',' + number(p.position.y) + ',' + number(p.position.z) + ',' + number(p.facing.x) +
             ',' + number(p.facing.y) + ',' + number(p.facing.z) + ',' + number(p.width) + ',' + number(p.height) + ',' +
             port_state_text(p.state) + ',' + traversal_text(p.traversal));

    line("");
    line("=== TUNNELS ===");
    line("id,map_id,start_port_id,end_port_id,code,name,maximum_slope_degrees,minimum_clearance_m,"
         "traversal_class,is_locked,generation_seed,tags,generation_json");
    for (const auto& wt : w.tunnels) {
        const Tunnel& t = wt.t;
        line(t.id + ',' + w.map_id + ',' + t.start_port_id + ',' + t.end_port_id + ',' + csv(t.code) + ',' +
             csv(t.name) + ',' + number(t.maximum_slope_degrees) + ',' + number(t.minimum_clearance) + ',' +
             traversal_text(t.traversal) + ',' + (t.locked ? "1" : "0") + ',' + std::to_string(t.generation.seed) +
             ',' + csv(tags_json(wt.tags)) + ',' + csv(wt.gen));
    }

    line("");
    line("=== TUNNEL_POINTS ===");
    line("tunnel_id,point_index,center_x,center_y,center_z,width_m,height_m");
    for (const auto& wt : w.tunnels)
        for (std::size_t i = 0; i < wt.t.centerline.size(); ++i) {
            const Profile& p = wt.t.centerline[i];
            line(wt.t.id + ',' + std::to_string(i) + ',' + number(p.center.x) + ',' + number(p.center.y) + ',' +
                 number(p.center.z) + ',' + number(p.width) + ',' + number(p.height));
        }
    return s;
}

Face positive_face(int axis) { return axis == 0 ? Face::positive_x : axis == 1 ? Face::positive_y : Face::positive_z; }
Face negative_face(int axis) { return axis == 0 ? Face::negative_x : axis == 1 ? Face::negative_y : Face::negative_z; }

void add_tunnel(World& w, const Tunnel& local, const Vec3& origin) {
    WorldTunnel t{world_tunnel(local, origin), {}, {}};
    const std::string what = label("tunnel", t.t.id, t.t.code);
    t.tags = canonical_tags(local.tags, what);
    t.gen = generation_json(local.generation, what);
    for (const auto& p : t.t.centerline) {
        const double hw = p.width / 2, hh = p.height / 2;
        w.include(p.center.x - hw, p.center.y - hw, p.center.z - hh);
        w.include(p.center.x + hw, p.center.y + hw, p.center.z + hh);
    }
    w.tunnels.push_back(std::move(t));
}

} // namespace

std::string to_moocow_csv(const std::vector<Region>& input, const std::string& map_id, const std::string& map_name) {
    World w;
    w.map_id = uuid(map_id, "map_id");
    check_required_text(map_name, "map_name");
    w.map_name = map_name;
    // maps.code is unique across the whole database (NOCASE); deriving it from
    // the fresh map ID keeps a new package from colliding with an existing map.
    w.map_code = "fbs-caves-" + w.map_id;
    if (input.empty()) fail(ErrorCode::invalid_input, "no regions were supplied.");

    // Canonical order: output must not depend on the order regions were
    // generated or passed in.
    std::vector<const Region*> regions;
    for (const auto& r : input) regions.push_back(&r);
    std::sort(regions.begin(), regions.end(), [](const Region* a, const Region* b) { return key_less(a->key, b->key); });

    std::size_t cavern_count = 0;
    std::set<std::string> region_ids;
    for (std::size_t i = 0; i < regions.size(); ++i) {
        const Region& r = *regions[i];
        const std::string where = "region " + key_text(r.key);
        if (r.version != schema_version)
            fail(ErrorCode::incompatible_version, where + " has schema version " + std::to_string(r.version) + ".");
        if (r.algorithm != algorithm_version)
            fail(ErrorCode::incompatible_version, where + " has algorithm '" + r.algorithm + "'.");
        finite(r.origin, where + " origin");
        finite(r.size, where + " size");
        if (!(r.size.x > 0 && r.size.y > 0 && r.size.z > 0)) fail(ErrorCode::invalid_input, where + " size must be positive.");
        if (i > 0 && !key_less(regions[i - 1]->key, r.key))
            fail(ErrorCode::invalid_input, where + " was supplied more than once.");
        if (!r.id.empty() && !region_ids.insert(r.id).second)
            fail(ErrorCode::invalid_input, where + " repeats region id '" + r.id + "'.");
        if (!same_settings(*regions.front(), r))
            fail(ErrorCode::incompatible_boundary,
                 where + " does not share world_seed, settings_id, size, version and algorithm with region " +
                     key_text(regions.front()->key) + ".");
        check_lattice(*regions.front(), r);
        validate(r);
        cavern_count += r.caverns.size();
    }
    if (cavern_count == 0) fail(ErrorCode::invalid_input, "the regions contain no caverns.");

    // Region records, offset to world coordinates.
    for (const Region* rp : regions) {
        const Region& r = *rp;
        const Vec3& o = r.origin;
        w.include(o.x, o.y, o.z);
        w.include(o.x + r.size.x, o.y + r.size.y, o.z + r.size.z);
        for (const auto& c : r.caverns) {
            WorldCavern wc{world_cavern(c, o), {}, {}};
            const std::string what = label("cavern", wc.c.id, wc.c.code);
            wc.tags = canonical_tags(c.tags, what);
            wc.gen = generation_json(c.generation, what);
            w.include(wc.c.bounds.minimum.x, wc.c.bounds.minimum.y, wc.c.bounds.minimum.z);
            w.include(wc.c.bounds.maximum.x, wc.c.bounds.maximum.y, wc.c.bounds.maximum.z);
            w.include(wc.c.center.x, wc.c.center.y, wc.c.center.z);
            for (const auto& p : wc.c.boundary) w.include_xy(p.x, p.y);
            w.caverns.push_back(std::move(wc));
        }
        for (const auto& f : r.floors) {
            WorldFloor wf{world_floor(f, o), {}};
            wf.tags = canonical_tags(f.tags, label("floor", wf.f.id, wf.f.code));
            if (!wf.f.boundary.empty()) {
                const Vec2 anchor = wf.f.boundary.front();
                for (const auto& p : wf.f.boundary) {
                    const double z = wf.f.base_z + wf.f.slope_x * (p.x - anchor.x) + wf.f.slope_y * (p.y - anchor.y);
                    if (!std::isfinite(z)) fail(ErrorCode::invalid_input, label("floor", wf.f.id, wf.f.code) + " surface is not finite.");
                    w.include(p.x, p.y, z - wf.f.variation_amplitude);
                    w.include_z(z + wf.f.variation_amplitude);
                }
            }
            w.include_z(wf.f.base_z);
            w.floors.push_back(std::move(wf));
        }
        for (const auto& c : r.cliffs) {
            w.cliffs.push_back(world_cliff(c, o));
            for (const auto& p : w.cliffs.back().boundary) w.include_xy(p.x, p.y);
        }
        for (const auto& p : r.ports) {
            w.ports.push_back(world_port(p, o));
            // Conservative opening extent around the port centre.
            const Port& q = w.ports.back();
            const double hw = q.width / 2, hh = q.height / 2;
            w.include(q.position.x - hw, q.position.y - hw, q.position.z - hh);
            w.include(q.position.x + hw, q.position.y + hw, q.position.z + hh);
        }
        for (const auto& t : r.tunnels) add_tunnel(w, t, o);
    }

    // Seams: every adjacent pair with a boundary on the shared face is joined
    // by the ordinary tunnel(s) stitch() returns. Every seam port on either
    // side must be an endpoint of a returned tunnel.
    std::map<std::tuple<std::int64_t, std::int64_t, std::int64_t>, const Region*> by_key;
    for (const Region* r : regions) by_key[{r->key.x, r->key.y, r->key.z}] = r;
    for (const Region* a : regions) {
        for (int axis = 0; axis < 3; ++axis) {
            RegionKey k = a->key;
            std::int64_t* component = axis == 0 ? &k.x : axis == 1 ? &k.y : &k.z;
            if (*component == std::numeric_limits<std::int64_t>::max()) continue;
            ++*component;
            const auto found = by_key.find({k.x, k.y, k.z});
            if (found == by_key.end()) continue;
            const Region* b = found->second;
            std::vector<std::string> seam_ports;
            for (const auto& boundary : a->boundaries)
                if (boundary.face == positive_face(axis)) seam_ports.push_back(uuid(boundary.port_id, "boundary port_id"));
            for (const auto& boundary : b->boundaries)
                if (boundary.face == negative_face(axis)) seam_ports.push_back(uuid(boundary.port_id, "boundary port_id"));
            if (seam_ports.empty()) continue;
            const std::vector<Tunnel> joined = stitch(*a, *b);
            std::unordered_set<std::string> endpoints;
            for (const auto& t : joined) {
                add_tunnel(w, t, a->origin);
                endpoints.insert(w.tunnels.back().t.start_port_id);
                endpoints.insert(w.tunnels.back().t.end_port_id);
            }
            for (const auto& port : seam_ports)
                if (!endpoints.count(port))
                    fail(ErrorCode::incompatible_boundary, "seam between regions " + key_text(a->key) + " and " +
                                                               key_text(b->key) + " leaves boundary port " + port +
                                                               " unconnected.");
        }
    }

    // Elevation bands: one per region Z layer, widened so together they cover
    // the whole map Z range without exclusive overlap.
    std::map<std::int64_t, std::pair<double, double>> layers;
    for (const Region* r : regions) layers.emplace(r->key.z, std::make_pair(r->origin.z, r->origin.z + r->size.z));
    if (layers.size() > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()))
        fail(ErrorCode::invalid_input, "too many elevation layers.");
    for (const auto& layer : layers) {
        Band b;
        b.id = uuid(stable_id("fbs.caves.moocow.elevation_band\n" + w.map_id + "\n" + std::to_string(layer.first)),
                    "stable_id elevation band");
        b.code = "fbs-z" + std::to_string(layer.first);
        b.name = "Region layer z=" + std::to_string(layer.first);
        b.min_z = layer.second.first;
        b.max_z = layer.second.second;
        b.order = static_cast<int>(w.bands.size());
        w.bands.push_back(b);
    }
    w.bands.front().min_z = w.min.z;
    w.bands.back().max_z = w.max.z;
    for (std::size_t i = 0; i + 1 < w.bands.size(); ++i) w.bands[i].max_z = w.bands[i + 1].min_z;
    for (const auto& b : w.bands)
        if (!(b.min_z <= b.max_z)) fail(ErrorCode::incompatible_boundary, "region Z layers overlap.");
    finite(w.min, "map minimum");
    finite(w.max, "map maximum");
    if (!std::isfinite(w.max.x - w.min.x) || !std::isfinite(w.max.y - w.min.y) || !std::isfinite(w.max.z - w.min.z))
        fail(ErrorCode::invalid_input, "map extent is not finite.");

    // Canonical record order is MOOCoW's storage order
    // (CaveDatasetRepository.LoadDatasetWithRevision): caverns and tunnels
    // ORDER BY id; floors, cliffs and ports per cavern in cavern order, each
    // ORDER BY id. IDs are TEXT with BINARY collation in the lower-case "D"
    // form. The CSV then equals MOOCoW's own export of the stored map, so its
    // revision content hash is reproducible from the package bytes.
    auto owned = [](const std::string& owner_a, const std::string& id_a, const std::string& owner_b,
                    const std::string& id_b) { return owner_a != owner_b ? owner_a < owner_b : id_a < id_b; };
    std::stable_sort(w.caverns.begin(), w.caverns.end(),
                     [](const WorldCavern& a, const WorldCavern& b) { return a.c.id < b.c.id; });
    std::stable_sort(w.floors.begin(), w.floors.end(), [&](const WorldFloor& a, const WorldFloor& b) {
        return owned(a.f.cavern_id, a.f.id, b.f.cavern_id, b.f.id);
    });
    std::stable_sort(w.cliffs.begin(), w.cliffs.end(),
                     [&](const Cliff& a, const Cliff& b) { return owned(a.cavern_id, a.id, b.cavern_id, b.id); });
    std::stable_sort(w.ports.begin(), w.ports.end(),
                     [&](const Port& a, const Port& b) { return owned(a.cavern_id, a.id, b.cavern_id, b.id); });
    std::stable_sort(w.tunnels.begin(), w.tunnels.end(),
                     [](const WorldTunnel& a, const WorldTunnel& b) { return a.t.id < b.t.id; });

    Issues issues;
    check_identity_and_codes(w, issues);
    check_records(w, issues);
    issues.raise();
    return emit(w);
}

} // namespace fbs::caves
