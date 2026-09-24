// Deterministic identities for FinalBuildCaves.
//
// stable_id(identity) returns a custom RFC 9562 version-8 UUID. The 128 payload
// bits come from Mixer128, a deterministic, platform-independent, domain-separated
// mixer (it is NOT cryptographic and must not be used for security):
//
//   state   a = 0x243f6a8885a308d3, b = 0x13198a2e03707344 (pi digits)
//   absorb  typed, length-prefixed 64-bit little-endian words; each word w does
//           a ^= w; b += w * 0x9e3779b97f4a7c15;
//           a = fmix(a + rotl(b, 17)); b = fmix(b ^ rotl(a, 29))
//           where fmix is the SplitMix64 finaliser (a bijection on 64 bits).
//   domain  every mixer first absorbs its domain string, so stable_id
//           ("fbs.caves.stable_id.v1"), seam derivation and RNG seeding never
//           share an input space.
//   finish  absorb the word count, then a = fmix(a ^ 0xa0761d6478bd642f);
//           b = fmix(b + a); a = fmix(a + rotl(b, 31)).
//   UUID    bytes = big-endian(a) || big-endian(b); byte 6 high nibble = 8
//           (version), byte 8 top bits = 10 (RFC variant). Lower-case 8-4-4-4-12.
//
// Identity strings are built from exact integers and exact (shortest round-trip)
// decimal forms of doubles. std::hash, coordinate rounding and process state are
// never used, so identities are stable across runs, threads and generation order.

#include "region_internal.hpp"

#include <charconv>
#include <cstring>
#include <limits>

namespace fbs::caves::detail {
namespace {
inline std::uint64_t rotl(std::uint64_t x, int r) { return (x << r) | (x >> (64 - r)); }
inline std::uint64_t fmix(std::uint64_t z) {
    z ^= z >> 30;
    z *= 0xbf58476d1ce4e5b9ULL;
    z ^= z >> 27;
    z *= 0x94d049bb133111ebULL;
    z ^= z >> 31;
    return z;
}
constexpr std::uint64_t tag_u64 = 0x7531;
constexpr std::uint64_t tag_str = 0x7332;
constexpr std::uint64_t tag_f64 = 0x6633;
} // namespace

Mixer128::Mixer128(const char* domain) : a_(0x243f6a8885a308d3ULL), b_(0x13198a2e03707344ULL) {
    str(std::string(domain));
}

void Mixer128::block(std::uint64_t w) {
    a_ ^= w;
    b_ += w * 0x9e3779b97f4a7c15ULL;
    a_ = fmix(a_ + rotl(b_, 17));
    b_ = fmix(b_ ^ rotl(a_, 29));
    ++count_;
}

Mixer128& Mixer128::bytes(const void* data, std::size_t size) {
    const auto* p = static_cast<const unsigned char*>(data);
    while (size > 0) {
        std::uint64_t w = 0;
        const std::size_t n = size < 8 ? size : 8;
        for (std::size_t i = 0; i < n; ++i) w |= static_cast<std::uint64_t>(p[i]) << (8 * i);
        block(w);
        p += n;
        size -= n;
    }
    return *this;
}

Mixer128& Mixer128::str(const std::string& text) {
    block(tag_str);
    block(static_cast<std::uint64_t>(text.size()));
    return bytes(text.data(), text.size());
}

Mixer128& Mixer128::u64(std::uint64_t value) {
    block(tag_u64);
    block(value);
    return *this;
}

Mixer128& Mixer128::f64(double value) {
    if (value == 0) value = 0; // unify -0 and +0
    std::uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof bits);
    block(tag_f64);
    block(bits);
    return *this;
}

Hash128 Mixer128::finish() const {
    Mixer128 copy = *this;
    copy.block(copy.count_);
    std::uint64_t a = fmix(copy.a_ ^ 0xa0761d6478bd642fULL);
    std::uint64_t b = fmix(copy.b_ + a);
    a = fmix(a + rotl(b, 31));
    return {a, b};
}

std::string format_uuid_v8(const Hash128& hash) {
    unsigned char bytes[16];
    for (int i = 0; i < 8; ++i) {
        bytes[i] = static_cast<unsigned char>(hash.hi >> (56 - 8 * i));
        bytes[8 + i] = static_cast<unsigned char>(hash.lo >> (56 - 8 * i));
    }
    bytes[6] = static_cast<unsigned char>((bytes[6] & 0x0f) | 0x80);
    bytes[8] = static_cast<unsigned char>((bytes[8] & 0x3f) | 0x80);
    static const char* hex = "0123456789abcdef";
    std::string out;
    out.reserve(36);
    for (int i = 0; i < 16; ++i) {
        if (i == 4 || i == 6 || i == 8 || i == 10) out.push_back('-');
        out.push_back(hex[bytes[i] >> 4]);
        out.push_back(hex[bytes[i] & 15]);
    }
    return out;
}

bool is_uuid(const std::string& text) {
    if (text.size() != 36) return false;
    bool nonzero = false;
    for (std::size_t i = 0; i < 36; ++i) {
        const char c = text[i];
        if (i == 8 || i == 13 || i == 18 || i == 23) {
            if (c != '-') return false;
            continue;
        }
        const bool digit = c >= '0' && c <= '9';
        const bool lower = c >= 'a' && c <= 'f';
        const bool upper = c >= 'A' && c <= 'F';
        if (!digit && !lower && !upper) return false;
        if (c != '0') nonzero = true;
    }
    return nonzero;
}

std::string uuid_key(const std::string& text) {
    std::string out = text;
    for (char& c : out)
        if (c >= 'A' && c <= 'Z') c = static_cast<char>(c - 'A' + 'a');
    return out;
}

std::uint64_t Rng::next() {
    std::uint64_t z = (state_ += 0x9e3779b97f4a7c15ULL);
    return fmix(z);
}

double Rng::uniform() { return static_cast<double>(next() >> 11) * 0x1.0p-53; }

std::uint32_t Rng::below(std::uint32_t n) {
    return static_cast<std::uint32_t>(((next() >> 32) * static_cast<std::uint64_t>(n)) >> 32);
}

std::string format_double(double value) {
    if (value == 0) value = 0;
    char buffer[64];
    const auto result = std::to_chars(buffer, buffer + sizeof buffer, value);
    return std::string(buffer, result.ptr);
}

bool parse_double(const std::string& text, double& out) {
    if (text.empty()) return false;
    const char* end = text.data() + text.size();
    const auto result = std::from_chars(text.data(), end, out);
    return result.ec == std::errc() && result.ptr == end && std::isfinite(out);
}

bool parse_u64(const std::string& text, std::uint64_t& out) {
    if (text.empty()) return false;
    const char* end = text.data() + text.size();
    const auto result = std::from_chars(text.data(), end, out);
    return result.ec == std::errc() && result.ptr == end;
}

std::string key_label(const RegionKey& key) {
    return std::to_string(key.x) + "_" + std::to_string(key.y) + "_" + std::to_string(key.z);
}

const char* face_label(Face face) {
    switch (face) {
    case Face::negative_x: return "nx";
    case Face::positive_x: return "px";
    case Face::negative_y: return "ny";
    case Face::positive_y: return "py";
    case Face::negative_z: return "nz";
    case Face::positive_z: return "pz";
    }
    return "??";
}

void Budget::check_cancelled() {
    if (options_.cancelled && options_.cancelled()) throw Error(ErrorCode::cancelled, "generation cancelled");
}

void Budget::charge(std::uint64_t units) {
    used_ = units > std::numeric_limits<std::uint64_t>::max() - used_ ? std::numeric_limits<std::uint64_t>::max()
                                                                       : used_ + units;
    if (used_ > options_.work_limit) throw Error(ErrorCode::budget, "regional work limit exceeded");
    since_check_ += units;
    if (since_check_ >= 4096) {
        since_check_ = 0;
        check_cancelled();
    }
}

bool is_generated(const Generation& generation) {
    const auto it = generation.parameters.find(generator_key);
    return it != generation.parameters.end() && it->second == algorithm_version;
}

} // namespace fbs::caves::detail

namespace fbs::caves {
std::string stable_id(const std::string& identity) {
    detail::Mixer128 mixer("fbs.caves.stable_id.v1");
    mixer.str(identity);
    return detail::format_uuid_v8(mixer.finish());
}
} // namespace fbs::caves
