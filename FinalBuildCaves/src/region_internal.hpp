#pragma once
// Internal declarations shared by the FinalBuildCaves regional generator.
// Nothing in this header is part of the public contract (include/fbs/caves.hpp).

#include "fbs/caves.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <utility>
#include <vector>

namespace fbs::caves::detail {

// ---------------------------------------------------------------------------
// Hard limits. Every loop and allocation in the regional generator is bounded
// by these constants or by the caller's Options. They are validated before any
// size-dependent allocation happens.
// ---------------------------------------------------------------------------
inline constexpr double max_region_extent = 100000.0;   // meters per axis
inline constexpr double min_region_extent = 8.0;        // meters per axis
inline constexpr double max_feature_size = 10000.0;     // meters (radii, widths, heights)
inline constexpr std::uint32_t max_chambers = 256;
inline constexpr std::uint32_t max_loops = 256;
inline constexpr std::uint32_t max_neighbors = 64;
inline constexpr std::size_t max_authored_records = 100000;   // per record kind
inline constexpr std::size_t max_polygon_vertices = 1024;
inline constexpr std::size_t max_total_vertices = 262144;     // all authored polygons and paths
inline constexpr std::size_t max_region_records = 1000000;    // any record kind in a Region passed to validate()
inline constexpr std::size_t max_tunnel_profiles = 100000;
inline constexpr std::size_t max_reserved_boxes = 10000;
inline constexpr std::size_t max_string_bytes = 4096;
inline constexpr std::size_t max_tags = 256;
inline constexpr std::size_t max_parameters = 1024;
inline constexpr std::uint64_t hard_max_support_points = 4000000;
inline constexpr std::uint64_t hard_max_support_edges = 64000000;
inline constexpr double geometry_tolerance = 1e-9;       // matches MOOCoW validator
inline constexpr double slope_margin_degrees = 0.25;     // generated edges keep this below the limit
inline constexpr double walkable_slope_degrees = 30.0;   // above this, generated passages are "difficult"
inline constexpr double pi = 3.14159265358979323846;

// Generated-record marker stored in Generation::parameters.
inline constexpr const char* generator_key = "fbs.generator";

// ---------------------------------------------------------------------------
// Deterministic 128-bit domain-separated mixer (non-cryptographic).
// ---------------------------------------------------------------------------
struct Hash128 {
    std::uint64_t hi = 0, lo = 0;
    bool operator==(const Hash128& o) const { return hi == o.hi && lo == o.lo; }
};

class Mixer128 {
public:
    explicit Mixer128(const char* domain);
    Mixer128& bytes(const void* data, std::size_t size);
    Mixer128& str(const std::string& text);     // length-prefixed
    Mixer128& u64(std::uint64_t value);
    Mixer128& i64(std::int64_t value) { return u64(static_cast<std::uint64_t>(value)); }
    Mixer128& f64(double value);                // exact bit pattern (+0 and -0 unified)
    Hash128 finish() const;

private:
    void block(std::uint64_t word);
    std::uint64_t a_, b_, count_ = 0;
};

std::string format_uuid_v8(const Hash128& hash);
// Accepts the canonical 8-4-4-4-12 hexadecimal form (either case), rejects nil.
bool is_uuid(const std::string& text);
// Lower-cased copy used for identity comparisons (MOOCoW Guid equality is case-insensitive).
std::string uuid_key(const std::string& text);

// SplitMix64 stream seeded from a Hash128. Portable: no std distributions.
class Rng {
public:
    explicit Rng(const Hash128& seed) : state_(seed.hi ^ (seed.lo * 0x9e3779b97f4a7c15ULL)) {}
    std::uint64_t next();
    double uniform();                              // [0, 1)
    double range(double lo, double hi) { return lo + (hi - lo) * uniform(); }
    std::uint32_t below(std::uint32_t n);          // [0, n), n > 0

private:
    std::uint64_t state_;
};

// ---------------------------------------------------------------------------
// Small vector helpers.
// ---------------------------------------------------------------------------
inline Vec3 add(Vec3 a, Vec3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline Vec3 sub(Vec3 a, Vec3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline Vec3 mul(Vec3 a, double s) { return {a.x * s, a.y * s, a.z * s}; }
inline double dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline double length(Vec3 a) { return std::sqrt(dot(a, a)); }
inline double distance(Vec3 a, Vec3 b) { return length(sub(a, b)); }
inline double horizontal(Vec3 a, Vec3 b) { return std::hypot(b.x - a.x, b.y - a.y); }
inline bool same(Vec3 a, Vec3 b) { return a.x == b.x && a.y == b.y && a.z == b.z; }
inline bool finite(Vec3 a) { return std::isfinite(a.x) && std::isfinite(a.y) && std::isfinite(a.z); }
inline bool finite(Vec2 a) { return std::isfinite(a.x) && std::isfinite(a.y); }

// MOOCoW stores port facings through UnitVector3, which divides by
// sqrt(x*x + y*y + z*z) (unfused, left to right). A facing survives storage
// unchanged only when it is a fixed point of that step. stable_unit iterates the
// identical step (products forced through memory, so no FMA contraction) until
// it reaches a fixed point; returns false when none is reached in 8 steps.
bool stable_unit(Vec3 direction, Vec3& out);
bool is_stable_unit(Vec3 v);

// MOOCoW slope rule: atan2(|dz|, horizontal) in degrees, 90 when horizontal <= tolerance.
double segment_slope_degrees(Vec3 a, Vec3 b);

// Point / segment / box distances. Boxes are axis aligned.
double point_box_distance(Vec3 p, const Bounds& box);
// Minimum distance between segment [a,b] and box (convex minimisation, conservative by 1e-6 m).
double segment_box_distance(Vec3 a, Vec3 b, const Bounds& box);
double segment_segment_distance(Vec3 p1, Vec3 q1, Vec3 p2, Vec3 q2);
double point_segment_distance(Vec3 p, Vec3 a, Vec3 b);
Bounds inflate(const Bounds& box, double margin);
bool boxes_overlap(const Bounds& a, const Bounds& b);

// Planar polygon helpers (implicitly closed, first vertex not repeated).
bool polygon_is_simple(const Polygon& polygon);          // MOOCoW rules: no consecutive duplicates, no crossings
bool polygon_contains(const Polygon& polygon, Vec2 p);   // boundary counts as inside (MOOCoW rule)
bool has_consecutive_duplicates(const std::vector<Vec2>& path);
double polygon_area(const Polygon& polygon);

// Capsule radius conservatively covering a width x height rectangular section.
inline double section_radius(double width, double height) { return 0.5 * std::hypot(width, height); }

// ---------------------------------------------------------------------------
// Exclusion volumes used while routing.
// ---------------------------------------------------------------------------
struct Exclusion {
    enum Kind { box, capsule } kind = box;
    Bounds bounds;               // box, or AABB of the capsule already including radius
    Vec3 a, b;                   // capsule segment
    double radius = 0;           // capsule radius
    double clearance = 0;        // extra clearance owned by this exclusion (authored reserved_clearance)
    std::string owner;           // cavern uuid_key for chamber boxes (legal endpoint access), else empty
};

// Distance from segment to the exclusion surface minus its own clearance.
double exclusion_distance(const Exclusion& e, Vec3 a, Vec3 b);

// ---------------------------------------------------------------------------
// Generation settings shared by adjacent regions.
// ---------------------------------------------------------------------------
struct Shared {
    std::uint64_t world_seed = 0, style_id = 0;
    Vec3 size;
    std::uint32_t chamber_count = 0, loop_count = 0, neighbors = 0;
    double chamber_radius_min = 0, chamber_radius_max = 0, chamber_height = 0;
    double tunnel_width = 0, tunnel_height = 0, maximum_slope_degrees = 0, sampling_radius = 0;
    bool shelves = true;
};

Shared shared_from_request(const Request& request);
// Canonical key/value strings (exact round-trip decimal via to_chars).
std::vector<std::pair<std::string, std::string>> shared_parameters(const Shared& shared);
// Parse back from Generation::parameters; returns false when keys are missing/malformed.
bool shared_from_parameters(const std::map<std::string, std::string>& parameters, Shared& out);
std::string settings_id(const Shared& shared);

std::string region_identity(const Shared& shared, const RegionKey& key);

// Derived dimensions common to generation, validation and stitching.
struct Derived {
    double gap = 0;              // minimum rock between generated passages / chambers
    double width_max = 0;        // max generated tunnel width
    double height_max = 0;       // max generated tunnel height
    double radius = 0;           // capsule radius of a generated tunnel
    double spacing = 0;          // support lattice spacing
    double gateway_radius = 0;   // gateway chamber nominal radius
    double gateway_height = 0;   // gateway chamber floor-to-ceiling
    double corridor = 0;         // horizontal gateway port to face distance (XY faces)
};
Derived derive(const Shared& shared);

// ---------------------------------------------------------------------------
// Canonical seams.
// ---------------------------------------------------------------------------
inline int face_axis(Face f) { return static_cast<int>(f) / 2; }
inline bool face_positive(Face f) { return (static_cast<int>(f) % 2) == 1; }
Face opposite(Face f);
// Neighbor key; returns false on signed 64-bit overflow.
bool neighbor_key(const RegionKey& key, Face face, RegionKey& out);

struct Seam {
    std::string key;           // shared UUID
    double width = 0, height = 0;
    Vec3 position;             // on the shared plane, in the local frame of the region asking
    Vec3 port;                 // gateway port centre for the region asking (local frame)
    Vec3 facing;               // outward horizontal facing of that gateway port
    double slope_degrees = 0;  // slope of the port->seam passage (0 for XY faces)
};
// Derives the seam of `face` for region `key`. Both neighbours obtain the same key,
// width, height and (after frame translation) position. Throws Error(constraint)
// when a Z face cannot honour the slope limit inside the region.
Seam derive_seam(const Shared& shared, const Derived& derived, const RegionKey& key, Face face);

// ---------------------------------------------------------------------------
// Misc.
// ---------------------------------------------------------------------------
std::string key_label(const RegionKey& key);          // "0_-1_2"
const char* face_label(Face face);                    // "nx", "px", ...
std::string format_double(double value);              // shortest round-trip (to_chars)
bool parse_double(const std::string& text, double& out);
bool parse_u64(const std::string& text, std::uint64_t& out);

// Cooperative cancellation / work accounting for the geometry phase.
class Budget {
public:
    explicit Budget(const Options& options) : options_(options) {}
    // Adds work units; throws Error(budget) past the limit, Error(cancelled) when requested.
    void charge(std::uint64_t units);
    void check_cancelled();
    std::uint64_t used() const { return used_; }
    std::uint64_t remaining() const { return used_ >= options_.work_limit ? 0 : options_.work_limit - used_; }
    const Options& options() const { return options_; }

private:
    const Options& options_;
    std::uint64_t used_ = 0;
    std::uint64_t since_check_ = 0;
};

// Input validation (region_validate.cpp). Throws Error(invalid_input / incompatible_version).
// The budget (when given) is charged for the work and polled for cancellation.
void validate_request(const Request& request, const Options& options, Budget* budget);
void check_region(const Region& region, Budget* budget);

// Is this record produced by the generator (marker parameter present)?
bool is_generated(const Generation& generation);

} // namespace fbs::caves::detail
