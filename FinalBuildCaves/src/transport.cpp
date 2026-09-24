#include "fbs/caves.hpp"
#include "transport_internal.hpp"
#include <nlohmann/json.hpp>
#include <cmath>
#include <limits>
#include <set>
#include <type_traits>

namespace fbs::caves {
namespace {
using Json = nlohmann::json;
constexpr std::size_t max_input_bytes = 16u * 1024u * 1024u;
void invalid(const std::string& message) { throw Error(ErrorCode::invalid_input, message); }
void keys(const Json& j, std::initializer_list<const char*> accepted) {
    if (!j.is_object()) invalid("Expected a JSON object");
    std::set<std::string> names;
    for (const auto* name : accepted) names.emplace(name);
    for (auto it = j.begin(); it != j.end(); ++it)
        if (!names.count(it.key())) invalid("Unknown field: " + it.key());
}
template<class T> T integer(const Json& j) {
    if (!j.is_number_integer()) invalid("Expected an integer");
    if (j.is_number_unsigned()) {
        const auto n = j.get<std::uint64_t>();
        if (n > static_cast<std::uint64_t>(std::numeric_limits<T>::max())) invalid("Integer out of range");
        return static_cast<T>(n);
    }
    const auto n = j.get<std::int64_t>();
    if constexpr (std::is_unsigned_v<T>) {
        if (n < 0 || static_cast<std::uint64_t>(n) > std::numeric_limits<T>::max()) invalid("Integer out of range");
    } else if (n < std::numeric_limits<T>::min() || n > std::numeric_limits<T>::max()) invalid("Integer out of range");
    return static_cast<T>(n);
}
double real(const Json& j) {
    if (!j.is_number()) invalid("Expected a finite number");
    const auto n = j.get<double>();
    if (!std::isfinite(n)) invalid("Expected a finite number");
    return n;
}
void finite_json(const Json& j) {
    if (j.is_number_float() && !std::isfinite(j.get<double>()))
        invalid("Cannot serialize a non-finite number");
    if (j.is_structured()) for (const auto& child : j) finite_json(child);
}
template<class T> T value(const Json& j) {
    if constexpr (std::is_same_v<T, double>) return real(j);
    else if constexpr (std::is_same_v<T, bool>) {
        if (!j.is_boolean()) invalid("Expected a boolean");
        return j.get<bool>();
    } else if constexpr (std::is_integral_v<T>) return integer<T>(j);
    else return j.get<T>();
}
template<class T> void read(const Json& j, const char* name, T& target) {
    const auto it = j.find(name);
    if (it != j.end()) target = value<T>(*it);
}
Json parse(const std::string& text) {
    if (text.size() > max_input_bytes) invalid("JSON request exceeds 16 MiB");
    // Bound nesting before the parser allocates a recursive object graph.
    return Json::parse(text, [seen=std::map<int,std::set<std::string>>{}](int depth, Json::parse_event_t event, Json& parsed) mutable {
        if (depth > 64) invalid("JSON nesting exceeds 64 levels");
        if (event == Json::parse_event_t::object_start) seen[depth+1].clear();
        if (event == Json::parse_event_t::key && !seen[depth].insert(parsed.get<std::string>()).second) invalid("Duplicate JSON field");
        return true;
    });
}
template<class T> const char* enum_name(T v, std::initializer_list<std::pair<T,const char*>> pairs) {
    for (const auto& p : pairs) if (p.first == v) return p.second;
    invalid("Unknown enum value"); return "";
}
template<class T> T enum_value(const Json& j, std::initializer_list<std::pair<T,const char*>> pairs) {
    if (!j.is_string()) invalid("Enum must be a string");
    for (const auto& p : pairs) if (j == p.second) return p.first;
    invalid("Unknown enum name"); return T{};
}
} // namespace

#define ENUM_CODEC(Type, ...) \
void to_json(Json& j, Type v) { j = enum_name<Type>(v, {__VA_ARGS__}); } \
void from_json(const Json& j, Type& v) { v = enum_value<Type>(j, {__VA_ARGS__}); }
ENUM_CODEC(Traversal, {Traversal::walkable,"walkable"}, {Traversal::difficult,"difficult"}, {Traversal::hazardous,"hazardous"}, {Traversal::restricted,"restricted"}, {Traversal::impassable,"impassable"})
ENUM_CODEC(PortState, {PortState::required,"required"}, {PortState::optional,"optional"}, {PortState::generated,"generated"}, {PortState::sealed,"sealed"})
ENUM_CODEC(Authoring, {Authoring::exact,"exact"}, {Authoring::constrained_procedural,"constrained_procedural"}, {Authoring::generated,"generated"})
ENUM_CODEC(Transition, {Transition::walkable,"walkable"}, {Transition::slope,"slope"}, {Transition::ramp,"ramp"}, {Transition::stairs,"stairs"}, {Transition::bridge,"bridge"}, {Transition::climbable_cliff,"climbable_cliff"}, {Transition::impassable_cliff,"impassable_cliff"}, {Transition::one_way_drop,"one_way_drop"}, {Transition::unconnected,"unconnected"})
ENUM_CODEC(Face, {Face::negative_x,"negative_x"}, {Face::positive_x,"positive_x"}, {Face::negative_y,"negative_y"}, {Face::positive_y,"positive_y"}, {Face::negative_z,"negative_z"}, {Face::positive_z,"positive_z"})
#undef ENUM_CODEC

void to_json(Json& j, const Vec2& v) { j = Json::array({v.x,v.y}); }
void from_json(const Json& j, Vec2& v) {
    if (!j.is_array() || j.size()!=2) invalid("Vec2 needs exactly two numbers");
    v = {real(j[0]),real(j[1])};
}
void to_json(Json& j, const Vec3& v) { j = Json::array({v.x,v.y,v.z}); }
void from_json(const Json& j, Vec3& v) {
    if (!j.is_array() || j.size()!=3) invalid("Vec3 needs exactly three numbers");
    v = {real(j[0]),real(j[1]),real(j[2])};
}
void to_json(Json& j, const RegionKey& v) { j = Json::array({v.x,v.y,v.z}); }
void from_json(const Json& j, RegionKey& v) {
    if (!j.is_array() || j.size()!=3) invalid("Region key needs exactly three integers");
    v = {integer<std::int64_t>(j[0]),integer<std::int64_t>(j[1]),integer<std::int64_t>(j[2])};
}
void to_json(Json& j, const Bounds& v) { j = Json{{"minimum",v.minimum},{"maximum",v.maximum}}; }
void from_json(const Json& j, Bounds& v) {
    keys(j,{"minimum","maximum"});
    if (!j.contains("minimum") || !j.contains("maximum")) invalid("Bounds need minimum and maximum");
    read(j,"minimum",v.minimum); read(j,"maximum",v.maximum);
}
void to_json(Json& j, const Generation& v) { j = Json{{"seed",v.seed},{"parameters",v.parameters}}; }
void from_json(const Json& j, Generation& v) { keys(j,{"seed","parameters"}); read(j,"seed",v.seed); read(j,"parameters",v.parameters); }

#define PUT(field) j[#field] = v.field
#define GET(field) read(j,#field,v.field)
void to_json(Json& j, const Cavern& v) {
    j=Json::object(); PUT(id);PUT(code);PUT(name);PUT(center);PUT(bounds);PUT(boundary);PUT(authoring);PUT(locked);PUT(reserved_clearance);PUT(generation);PUT(tags);
}
void from_json(const Json& j, Cavern& v) {
    keys(j,{"id","code","name","center","bounds","boundary","authoring","locked","reserved_clearance","generation","tags"});
    GET(id);GET(code);GET(name);GET(center);GET(bounds);GET(boundary);GET(authoring);GET(locked);GET(reserved_clearance);GET(generation);GET(tags);
}
void to_json(Json& j, const Floor& v) {
    j=Json::object();PUT(id);PUT(cavern_id);PUT(code);PUT(name);PUT(boundary);PUT(base_z);PUT(slope_x);PUT(slope_y);PUT(variation_amplitude);PUT(variation_seed);PUT(material);PUT(traversal);PUT(tags);
}
void from_json(const Json& j, Floor& v) {
    keys(j,{"id","cavern_id","code","name","boundary","base_z","slope_x","slope_y","variation_amplitude","variation_seed","material","traversal","tags"});
    GET(id);GET(cavern_id);GET(code);GET(name);GET(boundary);GET(base_z);GET(slope_x);GET(slope_y);GET(variation_amplitude);GET(variation_seed);GET(material);GET(traversal);GET(tags);
}
void to_json(Json& j, const Cliff& v) {
    j=Json::object();PUT(id);PUT(cavern_id);PUT(from_floor_id);PUT(to_floor_id);PUT(code);PUT(boundary);PUT(transition);PUT(height);
}
void from_json(const Json& j, Cliff& v) {
    keys(j,{"id","cavern_id","from_floor_id","to_floor_id","code","boundary","transition","height"});
    GET(id);GET(cavern_id);GET(from_floor_id);GET(to_floor_id);GET(code);GET(boundary);GET(transition);GET(height);
}
void to_json(Json& j, const Port& v) {
    j=Json::object();PUT(id);PUT(cavern_id);PUT(floor_id);PUT(code);PUT(name);PUT(position);PUT(facing);PUT(width);PUT(height);PUT(state);PUT(traversal);
}
void from_json(const Json& j, Port& v) {
    keys(j,{"id","cavern_id","floor_id","code","name","position","facing","width","height","state","traversal"});
    GET(id);GET(cavern_id);GET(floor_id);GET(code);GET(name);GET(position);GET(facing);GET(width);GET(height);GET(state);GET(traversal);
}
void to_json(Json& j, const Profile& v) { j=Json::object();PUT(center);PUT(width);PUT(height); }
void from_json(const Json& j, Profile& v) { keys(j,{"center","width","height"});GET(center);GET(width);GET(height); }
void to_json(Json& j, const Tunnel& v) {
    j=Json::object();PUT(id);PUT(start_port_id);PUT(end_port_id);PUT(code);PUT(name);PUT(centerline);PUT(maximum_slope_degrees);PUT(minimum_clearance);PUT(traversal);PUT(locked);PUT(generation);PUT(tags);PUT(source_node_ids);
}
void from_json(const Json& j, Tunnel& v) {
    keys(j,{"id","start_port_id","end_port_id","code","name","centerline","maximum_slope_degrees","minimum_clearance","traversal","locked","generation","tags","source_node_ids"});
    GET(id);GET(start_port_id);GET(end_port_id);GET(code);GET(name);GET(centerline);GET(maximum_slope_degrees);GET(minimum_clearance);GET(traversal);GET(locked);GET(generation);GET(tags);
    if (j.contains("source_node_ids")) {
        if (!j["source_node_ids"].is_array()) invalid("source_node_ids must be an array");
        for (const auto& n : j["source_node_ids"]) v.source_node_ids.push_back(integer<std::uint32_t>(n));
    }
}
void to_json(Json& j, const Boundary& v) { j=Json::object();PUT(key);PUT(port_id);PUT(face);PUT(position);PUT(width);PUT(height); }
void from_json(const Json& j, Boundary& v) { keys(j,{"key","port_id","face","position","width","height"});GET(key);GET(port_id);GET(face);GET(position);GET(width);GET(height); }
void to_json(Json& j, const Region& v) {
    j=Json::object();PUT(version);PUT(algorithm);PUT(id);PUT(settings_id);PUT(world_seed);PUT(key);PUT(size);PUT(origin);PUT(caverns);PUT(floors);PUT(cliffs);PUT(ports);PUT(tunnels);PUT(boundaries);PUT(work_used);PUT(support_points);
}
void from_json(const Json& j, Region& v) {
    keys(j,{"version","algorithm","id","settings_id","world_seed","key","size","origin","caverns","floors","cliffs","ports","tunnels","boundaries","work_used","support_points"});
    GET(version);GET(algorithm);GET(id);GET(settings_id);GET(world_seed);GET(key);GET(size);GET(origin);GET(caverns);GET(floors);GET(cliffs);GET(ports);GET(tunnels);GET(boundaries);GET(work_used);GET(support_points);
    if (v.version != schema_version || v.algorithm != algorithm_version) throw Error(ErrorCode::incompatible_version,"Unsupported region schema or algorithm");
}
void to_json(Json& j, const Request& v) {
    j=Json::object();PUT(version);PUT(world_seed);PUT(key);PUT(size);PUT(face_mask);PUT(chamber_count);PUT(loop_count);PUT(chamber_radius_min);PUT(chamber_radius_max);PUT(chamber_height);PUT(tunnel_width);PUT(tunnel_height);PUT(maximum_slope_degrees);PUT(sampling_radius);PUT(neighbors);PUT(shelves);PUT(style_id);PUT(reserved);PUT(authored_caverns);PUT(authored_floors);PUT(authored_cliffs);PUT(authored_ports);PUT(authored_tunnels);
}
void from_json(const Json& j, Request& v) {
    keys(j,{"version","world_seed","key","size","face_mask","chamber_count","loop_count","chamber_radius_min","chamber_radius_max","chamber_height","tunnel_width","tunnel_height","maximum_slope_degrees","sampling_radius","neighbors","shelves","style_id","reserved","authored_caverns","authored_floors","authored_cliffs","authored_ports","authored_tunnels"});
    GET(version);GET(world_seed);GET(key);GET(size);GET(face_mask);GET(chamber_count);GET(loop_count);GET(chamber_radius_min);GET(chamber_radius_max);GET(chamber_height);GET(tunnel_width);GET(tunnel_height);GET(maximum_slope_degrees);GET(sampling_radius);GET(neighbors);GET(shelves);GET(style_id);GET(reserved);GET(authored_caverns);GET(authored_floors);GET(authored_cliffs);GET(authored_ports);GET(authored_tunnels);
    if (v.version != schema_version) throw Error(ErrorCode::incompatible_version,"Unsupported request schema");
}
#undef PUT
#undef GET

std::string detail::serialize_validated_region(const Region& region) {
    try { const Json json(region);finite_json(json);return json.dump(2) + "\n"; }
    catch (const Json::exception& e) { invalid(e.what()); }
    throw Error(ErrorCode::internal,"Unreachable JSON path");
}
std::string to_json(const Region& region) { validate(region);return detail::serialize_validated_region(region); }
std::string request_to_json(const Request& request) {
    try { const Json json(request);finite_json(json);return json.dump(2) + "\n"; }
    catch (const Json::exception& e) { invalid(e.what()); }
    throw Error(ErrorCode::internal,"Unreachable JSON path");
}
Region region_from_json(const std::string& text) {
    try { auto region=parse(text).get<Region>();validate(region);return region; }
    catch (const Error&) { throw; }
    catch (const Json::exception& e) { invalid(e.what()); }
    throw Error(ErrorCode::internal,"Unreachable JSON path");
}
Request request_from_json(const std::string& text) {
    try { return parse(text).get<Request>(); }
    catch (const Error&) { throw; }
    catch (const Json::exception& e) { invalid(e.what()); }
    throw Error(ErrorCode::internal,"Unreachable JSON path");
}
} // namespace fbs::caves
