#pragma once
#include "caves_export.h"

#include <cstdint>
#include <functional>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace fbs::caves {
inline constexpr std::uint32_t schema_version = 1;
inline constexpr const char* algorithm_version = "fbs-caves-1";

struct Vec2 { double x = 0, y = 0; };
struct Vec3 { double x = 0, y = 0, z = 0; };
struct Bounds { Vec3 minimum, maximum; };
struct RegionKey { std::int64_t x = 0, y = 0, z = 0; };
using Polygon = std::vector<Vec2>;
enum class Face : std::uint32_t { negative_x, positive_x, negative_y, positive_y, negative_z, positive_z };
enum class Traversal { walkable, difficult, hazardous, restricted, impassable };
enum class PortState { required, optional, generated, sealed };
enum class Authoring { exact, constrained_procedural, generated };
enum class Transition { walkable, slope, ramp, stairs, bridge, climbable_cliff, impassable_cliff, one_way_drop, unconnected };
enum class ErrorCode { invalid_input, incompatible_version, cancelled, budget, no_route, constraint, incompatible_boundary, internal };
class FBS_CAVES_API Error : public std::runtime_error {
public:
    ErrorCode code;
    Error(ErrorCode code_, const std::string& message) : std::runtime_error(message), code(code_) {}
};
struct Generation {
    std::int32_t seed = 0;
    std::map<std::string, std::string> parameters;
};
struct Cavern {
    std::string id, code, name;
    Vec3 center;
    Bounds bounds;
    Polygon boundary;
    Authoring authoring = Authoring::generated;
    bool locked = false;
    double reserved_clearance = 0;
    Generation generation;
    std::vector<std::string> tags;
};
struct Floor {
    std::string id, cavern_id, code, name;
    Polygon boundary;
    double base_z = 0, slope_x = 0, slope_y = 0, variation_amplitude = 0;
    std::int32_t variation_seed = 0;
    std::string material = "limestone";
    Traversal traversal = Traversal::walkable;
    std::vector<std::string> tags;
};
struct Cliff {
    std::string id, cavern_id, from_floor_id, to_floor_id, code;
    std::vector<Vec2> boundary;
    Transition transition = Transition::unconnected;
    double height = 0;
};
struct Port {
    std::string id, cavern_id, floor_id, code, name;
    Vec3 position, facing;
    double width = 4, height = 3;
    PortState state = PortState::generated;
    Traversal traversal = Traversal::walkable;
};
struct Profile { Vec3 center; double width = 4, height = 3; };
struct Tunnel {
    std::string id, start_port_id, end_port_id, code, name;
    std::vector<Profile> centerline;
    double maximum_slope_degrees = 60, minimum_clearance = 3;
    Traversal traversal = Traversal::difficult;
    bool locked = false;
    Generation generation;
    std::vector<std::string> tags;
    // Native skeleton identities are provenance, never persistent MOOCoW columns.
    std::vector<std::uint32_t> source_node_ids;
};
struct Boundary {
    // One common key/position/dimensions at the shared face; each side has its own
    // existing cavern port. stitch() returns the connecting ordinary tunnel.
    std::string key, port_id;
    Face face = Face::negative_x;
    Vec3 position;
    double width = 4, height = 3;
};
struct Region {
    std::uint32_t version = schema_version;
    std::string algorithm = algorithm_version, id, settings_id;
    std::uint64_t world_seed = 0;
    RegionKey key;
    Vec3 size{512, 512, 128}, origin;
    std::vector<Cavern> caverns;
    std::vector<Floor> floors;
    std::vector<Cliff> cliffs;
    std::vector<Port> ports;
    std::vector<Tunnel> tunnels;
    std::vector<Boundary> boundaries;
    std::uint64_t work_used = 0, support_points = 0;
};
struct Request {
    std::uint32_t version = schema_version;
    std::uint64_t world_seed = 1;
    RegionKey key;
    Vec3 size{512, 512, 128};
    // Cardinal horizontal faces by default; explicit Z-face support may reject
    // geometrically infeasible slope requests instead of fabricating walkability.
    std::uint32_t face_mask = 15;
    std::uint32_t chamber_count = 8, loop_count = 1;
    double chamber_radius_min = 8, chamber_radius_max = 18;
    double chamber_height = 12, tunnel_width = 4, tunnel_height = 3;
    double maximum_slope_degrees = 60;
    double sampling_radius = 0.08;
    std::uint32_t neighbors = 20;
    bool shelves = true;
    // World-generation configuration identity must be shared by adjacent jobs.
    std::uint64_t style_id = 1;
    std::vector<Bounds> reserved;
    // Accepted authored records are copied without mutation. All coordinates
    // are local meters, +Z up. IDs must be nonzero UUID strings for the adapter.
    std::vector<Cavern> authored_caverns;
    std::vector<Floor> authored_floors;
    std::vector<Cliff> authored_cliffs;
    std::vector<Port> authored_ports;
    std::vector<Tunnel> authored_tunnels;
};
struct Options {
    std::uint64_t work_limit = 100000000;
    std::uint64_t max_points = 50000, max_edges = 2000000;
    std::function<bool()> cancelled;
};

// Strong failure guarantee: a return value exists only after validation. No
// application files, stdout, database or shared world state are touched.
FBS_CAVES_API Region generate(const Request&, const Options& = {});
FBS_CAVES_API void validate(const Region&);
// Inputs may be supplied in either order. Output coordinates are local to the
// first region. Only a matching, adjacent face is accepted; unrelated IDs stay.
FBS_CAVES_API std::vector<Tunnel> stitch(const Region& first, const Region& second);
// Deterministic custom UUID (version 8) from a domain-separated identity string.
FBS_CAVES_API std::string stable_id(const std::string& identity);

// Transport and existing-format adapter. Authored record fields survive the
// JSON transport; manifest metadata never requires MOOCoW schema changes.
FBS_CAVES_API std::string to_json(const Region&);
FBS_CAVES_API Region region_from_json(const std::string&);
FBS_CAVES_API Request request_from_json(const std::string&);
FBS_CAVES_API std::string request_to_json(const Request&);
FBS_CAVES_API std::string to_moocow_csv(const std::vector<Region>& regions,
                          const std::string& map_id,
                          const std::string& map_name = "Generated caves");
} // namespace fbs::caves
