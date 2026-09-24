// TEST-ONLY stand-ins for validate(), stitch() and stable_id() so the MOOCoW
// adapter can be compiled and exercised on its own before (or without) the
// regional library. Never link this file together with src/region*.cpp; the
// CMake adapter test links the real fbs::caves library instead.
//
// Build: see FinalBuildCaves/tests/adapter_acceptance.py (--standalone).

#include "fbs/caves.hpp"

#include <cstdio>
#include <string>

namespace fbs::caves {

void validate(const Region& region) {
    if (region.version != schema_version) throw Error(ErrorCode::incompatible_version, "stub validate: version");
}

std::string stable_id(const std::string& identity) {
    // FNV-1a based, two lanes; formatted as a version-8 / RFC 9562 variant UUID.
    std::uint64_t h1 = 1469598103934665603ull, h2 = 1099511628211ull ^ 0x9e3779b97f4a7c15ull;
    for (unsigned char c : identity) {
        h1 = (h1 ^ c) * 1099511628211ull;
        h2 = (h2 ^ c) * 0x100000001b3ull + 0x632be59bd9b4e019ull;
    }
    h1 = (h1 & 0xffffffffffff0fffull) | 0x0000000000008000ull;
    h2 = (h2 & 0x3fffffffffffffffull) | 0x8000000000000000ull;
    char text[37];
    std::snprintf(text, sizeof text, "%08llx-%04llx-%04llx-%04llx-%012llx",
                  static_cast<unsigned long long>(h1 >> 32), static_cast<unsigned long long>((h1 >> 16) & 0xffff),
                  static_cast<unsigned long long>(h1 & 0xffff), static_cast<unsigned long long>(h2 >> 48),
                  static_cast<unsigned long long>(h2 & 0xffffffffffffull));
    return text;
}

namespace {
Face opposite(Face f) {
    switch (f) {
    case Face::negative_x: return Face::positive_x;
    case Face::positive_x: return Face::negative_x;
    case Face::negative_y: return Face::positive_y;
    case Face::positive_y: return Face::negative_y;
    case Face::negative_z: return Face::positive_z;
    case Face::positive_z: return Face::negative_z;
    }
    return f;
}

const Port* find_port(const Region& r, const std::string& id) {
    for (const auto& p : r.ports)
        if (p.id == id) return &p;
    return nullptr;
}
} // namespace

// Joins boundaries with equal keys on facing sides: first port -> shared face
// point -> second port, in first-region-local coordinates.
std::vector<Tunnel> stitch(const Region& first, const Region& second) {
    if (first.world_seed != second.world_seed || first.settings_id != second.settings_id)
        throw Error(ErrorCode::incompatible_boundary, "stub stitch: settings differ");
    const Vec3 shift{second.origin.x - first.origin.x, second.origin.y - first.origin.y,
                     second.origin.z - first.origin.z};
    std::vector<Tunnel> out;
    for (const auto& a : first.boundaries) {
        const Boundary* match = nullptr;
        for (const auto& b : second.boundaries)
            if (b.key == a.key && b.face == opposite(a.face)) match = &b;
        const Port* pa = find_port(first, a.port_id);
        const Port* pb = match ? find_port(second, match->port_id) : nullptr;
        if (!match || !pa || !pb) {
            // Only a boundary facing the other region needs a partner.
            const bool faces_second =
                (a.face == Face::positive_x && second.key.x == first.key.x + 1) ||
                (a.face == Face::negative_x && second.key.x == first.key.x - 1) ||
                (a.face == Face::positive_y && second.key.y == first.key.y + 1) ||
                (a.face == Face::negative_y && second.key.y == first.key.y - 1) ||
                (a.face == Face::positive_z && second.key.z == first.key.z + 1) ||
                (a.face == Face::negative_z && second.key.z == first.key.z - 1);
            if (faces_second) throw Error(ErrorCode::incompatible_boundary, "stub stitch: unmatched seam " + a.key);
            continue;
        }
        if (a.width != match->width || a.height != match->height)
            throw Error(ErrorCode::incompatible_boundary, "stub stitch: seam dimensions differ " + a.key);
        Tunnel t;
        t.id = stable_id("stub-seam\n" + a.key);
        t.start_port_id = pa->id;
        t.end_port_id = pb->id;
        t.code = "seam-" + a.key;
        t.name = "Seam " + a.key;
        t.centerline = {{pa->position, a.width, a.height},
                        {a.position, a.width, a.height},
                        {{pb->position.x + shift.x, pb->position.y + shift.y, pb->position.z + shift.z},
                         a.width, a.height}};
        t.maximum_slope_degrees = 45;
        t.minimum_clearance = a.height;
        t.traversal = Traversal::walkable;
        t.generation.parameters["seam"] = a.key;
        out.push_back(t);
    }
    // Unmatched boundaries on the second region's facing side are inconsistent too.
    for (const auto& b : second.boundaries) {
        bool matched = false;
        for (const auto& a : first.boundaries) matched = matched || (a.key == b.key && a.face == opposite(b.face));
        const bool faces_first =
            (b.face == Face::negative_x && second.key.x == first.key.x + 1) ||
            (b.face == Face::positive_x && second.key.x == first.key.x - 1) ||
            (b.face == Face::negative_y && second.key.y == first.key.y + 1) ||
            (b.face == Face::positive_y && second.key.y == first.key.y - 1) ||
            (b.face == Face::negative_z && second.key.z == first.key.z + 1) ||
            (b.face == Face::positive_z && second.key.z == first.key.z - 1);
        if (faces_first && !matched) throw Error(ErrorCode::incompatible_boundary, "stub stitch: unmatched seam " + b.key);
    }
    return out;
}

} // namespace fbs::caves
