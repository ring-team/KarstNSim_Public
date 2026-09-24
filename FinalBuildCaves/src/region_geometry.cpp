// Geometry helpers for the regional generator. All functions are pure and
// bounded (no allocation proportional to anything but their inputs).

#include "region_internal.hpp"

#include <algorithm>

namespace fbs::caves::detail {

namespace {
double mul_unfused(double a, double b) {
    volatile double r = a * b;
    return r;
}
} // namespace

bool stable_unit(Vec3 v, Vec3& out) {
    for (int i = 0; i < 8; ++i) {
        const double m = std::sqrt(mul_unfused(v.x, v.x) + mul_unfused(v.y, v.y) + mul_unfused(v.z, v.z));
        if (!(m > 1e-12) || !std::isfinite(m)) return false;
        const Vec3 n{v.x / m, v.y / m, v.z / m};
        if (same(n, v)) {
            out = v;
            return true;
        }
        v = n;
    }
    return false;
}

bool is_stable_unit(Vec3 v) {
    Vec3 out;
    return stable_unit(v, out) && same(out, v);
}

double segment_slope_degrees(Vec3 a, Vec3 b) {
    const double h = horizontal(a, b);
    if (h <= geometry_tolerance) return 90.0;
    return std::atan2(std::fabs(b.z - a.z), h) * 180.0 / pi;
}

double point_box_distance(Vec3 p, const Bounds& box) {
    const double dx = std::max({box.minimum.x - p.x, 0.0, p.x - box.maximum.x});
    const double dy = std::max({box.minimum.y - p.y, 0.0, p.y - box.maximum.y});
    const double dz = std::max({box.minimum.z - p.z, 0.0, p.z - box.maximum.z});
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

double segment_box_distance(Vec3 a, Vec3 b, const Bounds& box) {
    // The distance from a moving point to a convex set is convex along the
    // segment, so golden-section search converges to the global minimum.
    const Vec3 d = sub(b, a);
    auto f = [&](double t) { return point_box_distance(add(a, mul(d, t)), box); };
    double lo = 0, hi = 1;
    const double g = 0.6180339887498949;
    double x1 = hi - g * (hi - lo), x2 = lo + g * (hi - lo);
    double f1 = f(x1), f2 = f(x2);
    for (int i = 0; i < 80; ++i) {
        if (f1 <= f2) {
            hi = x2; x2 = x1; f2 = f1;
            x1 = hi - g * (hi - lo); f1 = f(x1);
        } else {
            lo = x1; x1 = x2; f1 = f2;
            x2 = lo + g * (hi - lo); f2 = f(x2);
        }
    }
    const double best = std::min({f(0), f(1), f1, f2});
    // Conservative by a tiny margin covering the bracket width.
    return std::max(0.0, best - 1e-6 * (1.0 + length(d)));
}

double point_segment_distance(Vec3 p, Vec3 a, Vec3 b) {
    const Vec3 d = sub(b, a);
    const double dd = dot(d, d);
    double t = dd > 0 ? dot(sub(p, a), d) / dd : 0;
    t = std::clamp(t, 0.0, 1.0);
    return distance(p, add(a, mul(d, t)));
}

double segment_segment_distance(Vec3 p1, Vec3 q1, Vec3 p2, Vec3 q2) {
    // Ericson, Real-Time Collision Detection, 5.1.9.
    const Vec3 d1 = sub(q1, p1), d2 = sub(q2, p2), r = sub(p1, p2);
    const double a = dot(d1, d1), e = dot(d2, d2), f = dot(d2, r);
    const double eps = 1e-18;
    double s = 0, t = 0;
    if (a <= eps && e <= eps) return distance(p1, p2);
    if (a <= eps) {
        t = std::clamp(f / e, 0.0, 1.0);
    } else {
        const double c = dot(d1, r);
        if (e <= eps) {
            s = std::clamp(-c / a, 0.0, 1.0);
        } else {
            const double b = dot(d1, d2);
            const double denom = a * e - b * b;
            s = denom > eps ? std::clamp((b * f - c * e) / denom, 0.0, 1.0) : 0.0;
            t = (b * s + f) / e;
            if (t < 0) {
                t = 0;
                s = std::clamp(-c / a, 0.0, 1.0);
            } else if (t > 1) {
                t = 1;
                s = std::clamp((b - c) / a, 0.0, 1.0);
            }
        }
    }
    const double exact = distance(add(p1, mul(d1, s)), add(p2, mul(d2, t)));
    // Guard against cancellation in near-parallel cases with endpoint checks.
    return std::min({exact, point_segment_distance(p1, p2, q2), point_segment_distance(q1, p2, q2),
                     point_segment_distance(p2, p1, q1), point_segment_distance(q2, p1, q1)});
}

Bounds inflate(const Bounds& box, double m) {
    return {{box.minimum.x - m, box.minimum.y - m, box.minimum.z - m},
            {box.maximum.x + m, box.maximum.y + m, box.maximum.z + m}};
}

bool boxes_overlap(const Bounds& a, const Bounds& b) {
    return a.minimum.x <= b.maximum.x && a.maximum.x >= b.minimum.x && a.minimum.y <= b.maximum.y &&
           a.maximum.y >= b.minimum.y && a.minimum.z <= b.maximum.z && a.maximum.z >= b.minimum.z;
}

namespace {
double orientation(Vec2 a, Vec2 b, Vec2 c) { return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x); }
bool bbox_overlap(Vec2 a, Vec2 b, Vec2 c, Vec2 d) {
    const double t = geometry_tolerance;
    return std::max(std::min(a.x, b.x), std::min(c.x, d.x)) <= std::min(std::max(a.x, b.x), std::max(c.x, d.x)) + t &&
           std::max(std::min(a.y, b.y), std::min(c.y, d.y)) <= std::min(std::max(a.y, b.y), std::max(c.y, d.y)) + t;
}
// Same predicate as the MOOCoW validator.
bool segments_intersect(Vec2 a, Vec2 b, Vec2 c, Vec2 d) {
    return orientation(a, b, c) * orientation(a, b, d) <= geometry_tolerance &&
           orientation(c, d, a) * orientation(c, d, b) <= geometry_tolerance && bbox_overlap(a, b, c, d);
}
bool on_segment(Vec2 a, Vec2 b, Vec2 p) {
    return std::fabs(orientation(a, b, p)) <= geometry_tolerance && bbox_overlap(a, b, p, p);
}
bool same2(Vec2 a, Vec2 b) { return a.x == b.x && a.y == b.y; }
} // namespace

bool has_consecutive_duplicates(const std::vector<Vec2>& path) {
    for (std::size_t i = 1; i < path.size(); ++i)
        if (same2(path[i - 1], path[i])) return true;
    return false;
}

bool polygon_is_simple(const Polygon& p) {
    const std::size_t n = p.size();
    if (n < 3) return false;
    if (has_consecutive_duplicates(p) || same2(p[0], p[n - 1])) return false;
    std::size_t distinct = 0;
    for (std::size_t i = 0; i < n && distinct < 3; ++i) {
        bool seen = false;
        for (std::size_t j = 0; j < i; ++j)
            if (same2(p[i], p[j])) { seen = true; break; }
        if (!seen) ++distinct;
    }
    if (distinct < 3) return false;
    for (std::size_t i = 0; i < n; ++i) {
        const Vec2 a1 = p[i], a2 = p[(i + 1) % n];
        for (std::size_t j = i + 1; j < n; ++j) {
            if (j == i + 1 || (i == 0 && j == n - 1)) continue;
            if (segments_intersect(a1, a2, p[j], p[(j + 1) % n])) return false;
        }
    }
    return std::fabs(polygon_area(p)) > geometry_tolerance;
}

bool polygon_contains(const Polygon& polygon, Vec2 point) {
    bool inside = false;
    const std::size_t n = polygon.size();
    for (std::size_t i = 0, j = n - 1; i < n; j = i++) {
        const Vec2 a = polygon[i], b = polygon[j];
        if (on_segment(a, b, point)) return true;
        if ((a.y > point.y) != (b.y > point.y) && point.x < (b.x - a.x) * (point.y - a.y) / (b.y - a.y) + a.x)
            inside = !inside;
    }
    return inside;
}

double polygon_area(const Polygon& p) {
    double s = 0;
    for (std::size_t i = 0, n = p.size(); i < n; ++i) {
        const Vec2 a = p[i], b = p[(i + 1) % n];
        s += a.x * b.y - b.x * a.y;
    }
    return 0.5 * s;
}

double exclusion_distance(const Exclusion& e, Vec3 a, Vec3 b) {
    if (e.kind == Exclusion::box) return segment_box_distance(a, b, e.bounds) - e.clearance;
    return segment_segment_distance(a, b, e.a, e.b) - e.radius - e.clearance;
}

} // namespace fbs::caves::detail
