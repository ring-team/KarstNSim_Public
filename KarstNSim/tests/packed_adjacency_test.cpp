// Standalone tests for KarstNSim::PackedAdjacency. No external test library.
// Checks stay active in Release builds (no reliance on assert).
//
// Build example:
//   g++ -std=c++17 -Wall -Wextra -Wpedantic -I KarstNSim/include KarstNSim/tests/packed_adjacency_test.cpp -o packed_adjacency_test

#include "KarstNSim/packed_adjacency.h"

#include <cstdint>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

using KarstNSim::PackedAdjacency;

namespace {

	int g_failures = 0;

#define CHECK(cond) \
	do { \
		if (!(cond)) { \
			std::fprintf(stderr, "%s:%d: CHECK failed: %s\n", __FILE__, __LINE__, #cond); \
			++g_failures; \
		} \
	} while (0)

#define CHECK_THROWS(expr, ExType) \
	do { \
		bool thrown_ = false; \
		try { expr; } \
		catch (const ExType&) { thrown_ = true; } \
		catch (...) {} \
		if (!thrown_) { \
			std::fprintf(stderr, "%s:%d: expected %s from: %s\n", __FILE__, __LINE__, #ExType, #expr); \
			++g_failures; \
		} \
	} while (0)

	// Compile-time constness checks.
	using CEdge = decltype(std::declval<const PackedAdjacency&>()(0, 0));
	using MEdge = decltype(std::declval<PackedAdjacency&>()(0, 0));
	static_assert(std::is_same<CEdge, PackedAdjacency::ConstEdgeRef>::value, "const operator() yields ConstEdgeRef");
	static_assert(std::is_same<MEdge, PackedAdjacency::EdgeRef>::value, "mutable operator() yields EdgeRef");
	static_assert(std::is_same<decltype(std::declval<const PackedAdjacency&>()[0][0]), PackedAdjacency::ConstEdgeRef>::value, "const row proxy");
	static_assert(std::is_same<decltype(std::declval<const PackedAdjacency&>().row(0)[0]), PackedAdjacency::ConstEdgeRef>::value, "const row()");
	static_assert(std::is_same<decltype(std::declval<const PackedAdjacency&>().at(0, 0)), PackedAdjacency::ConstEdgeRef>::value, "const at()");
	static_assert(std::is_same<decltype(std::declval<PackedAdjacency&>()[0][0]), PackedAdjacency::EdgeRef>::value, "mutable row proxy");
	static_assert(std::is_same<decltype(std::declval<CEdge>().target), const std::int32_t&>::value, "const target ref");
	static_assert(std::is_same<decltype(std::declval<MEdge>().target), std::int32_t&>::value, "mutable target ref");
	static_assert(std::is_same<decltype(std::declval<CEdge>().weight[0]), const float&>::value, "const weight element");
	static_assert(std::is_same<decltype(std::declval<MEdge>().weight[0]), float&>::value, "mutable weight element");
	static_assert(std::is_same<decltype(std::declval<PackedAdjacency::ConstWeightView>().begin()), const float*>::value, "const weight iter");
	static_assert(std::is_nothrow_move_constructible<PackedAdjacency>::value, "nothrow move");
	static_assert(std::is_nothrow_move_assignable<PackedAdjacency>::value, "nothrow move assign");
	static_assert(std::is_copy_constructible<PackedAdjacency>::value, "copyable");
	static_assert(std::is_trivially_copyable<PackedAdjacency::WeightView>::value, "lightweight view");
	static_assert(std::is_trivially_copyable<PackedAdjacency::ConstWeightView>::value, "lightweight view");

	void test_empty() {
		PackedAdjacency g;
		CHECK(g.size() == 0);
		CHECK(g.cols() == 0);
		CHECK(g.channels() == 0);
		CHECK(g.storage_bytes() == 0);
		CHECK_THROWS(g.at(0, 0), std::out_of_range);
		CHECK_THROWS(g.set(0, 0, -1, {}), std::out_of_range);
	}

	void test_multichannel_read_write() {
		PackedAdjacency g;
		g.reset(4, 3, 2);
		CHECK(g.size() == 4 && g.cols() == 3 && g.channels() == 2);

		// Fresh slots: unused, zero costs.
		for (std::size_t r = 0; r < g.size(); ++r) {
			CHECK(g.row(r).size() == 3);
			for (std::size_t s = 0; s < g.cols(); ++s) {
				auto e = g(r, s);
				CHECK(e.target == -1);
				CHECK(e.weight.size() == 2);
				for (float w : e.weight) CHECK(w == 0.0f);
			}
		}

		g.set(1, 2, 3, { 1.5f, 2.5f });
		CHECK(g(1, 2).target == 3);
		CHECK(g(1, 2).weight[0] == 1.5f && g(1, 2).weight[1] == 2.5f);
		CHECK(g[1][2].target == 3);
		CHECK(g.row(1)[2].weight[1] == 2.5f);
		CHECK(g.at(1, 2).weight.at(0) == 1.5f);
		CHECK_THROWS(g.at(1, 2).weight.at(2), std::out_of_range);

		// Neighbouring slots untouched (edge-major layout sanity).
		CHECK(g(1, 1).target == -1 && g(1, 1).weight[1] == 0.0f);
		CHECK(g(2, 0).target == -1 && g(2, 0).weight[0] == 0.0f);

		// Mutate through references.
		auto e = g[0][0];
		e.target = 2;
		e.weight[1] = 7.0f;
		for (float& w : g(3, 1).weight) w = 9.0f;
		CHECK(g(0, 0).target == 2 && g(0, 0).weight[0] == 0.0f && g(0, 0).weight[1] == 7.0f);
		CHECK(g(3, 1).weight[0] == 9.0f && g(3, 1).weight[1] == 9.0f);

		g.set_weights(1, 2, { -1.0f, -2.0f });
		CHECK(g(1, 2).target == 3);
		CHECK(g(1, 2).weight[0] == -1.0f && g(1, 2).weight[1] == -2.0f);

		// Unset a slot.
		g.set(1, 2, -1, { 0.0f, 0.0f });
		CHECK(g(1, 2).target == -1);

		// Const access.
		const PackedAdjacency& cg = g;
		CHECK(cg(0, 0).target == 2);
		CHECK(cg[0][0].weight[1] == 7.0f);
		CHECK(cg.row(3).size() == 3);
		float sum = 0.0f;
		for (float w : cg.at(3, 1).weight) sum += w;
		CHECK(sum == 18.0f);
		PackedAdjacency::ConstWeightView cv = g(3, 1).weight; // mutable-to-const conversion
		CHECK(cv.size() == 2 && cv[0] == 9.0f);
	}

	void test_zero_channels() {
		PackedAdjacency g;
		g.reset(3, 2, 0);
		CHECK(g.channels() == 0);
		auto e = g(2, 1);
		CHECK(e.weight.size() == 0 && e.weight.empty());
		CHECK(e.weight.data() == nullptr);
		CHECK(e.weight.begin() == e.weight.end());
		int iterations = 0;
		for (float w : e.weight) { (void)w; ++iterations; }
		CHECK(iterations == 0);
		const PackedAdjacency& cg = g;
		CHECK(cg(0, 0).weight.begin() == cg(0, 0).weight.end());
		g.set(0, 1, 2, {});
		CHECK(g(0, 1).target == 2);
		g.set_weights(0, 1, {});
		CHECK_THROWS(g.set(0, 1, 2, { 1.0f }), std::invalid_argument);
		CHECK(g.storage_bytes() == 3 * 2 * sizeof(std::int32_t));

		// Zero columns / zero rows.
		g.reset(5, 0, 2);
		CHECK(g.size() == 5 && g.cols() == 0 && g.row(4).size() == 0);
		CHECK_THROWS(g.at(0, 0), std::out_of_range);
		g.reset(0, 4, 2);
		CHECK(g.size() == 0 && g.cols() == 4);
		CHECK(g.storage_bytes() == 0);
	}

	void test_reset_shape() {
		PackedAdjacency g;
		g.reset(2, 2, 3);
		g.set(0, 0, 1, { 1.0f, 2.0f, 3.0f });
		g.reset(5, 4, 1);
		CHECK(g.size() == 5 && g.cols() == 4 && g.channels() == 1);
		for (std::size_t r = 0; r < 5; ++r)
			for (std::size_t s = 0; s < 4; ++s)
				CHECK(g(r, s).target == -1 && g(r, s).weight.size() == 1 && g(r, s).weight[0] == 0.0f);
		g.set(4, 3, 0, { 5.0f });
		CHECK(g(4, 3).target == 0 && g(4, 3).weight[0] == 5.0f);

		// Shrinking releases storage.
		g.reset(1, 1, 1);
		CHECK(g.storage_bytes() == sizeof(std::int32_t) + sizeof(float));
		CHECK(g(0, 0).target == -1);
	}

	void test_bad_inputs() {
		PackedAdjacency g;
		g.reset(3, 2, 2);
		g.set(0, 0, 1, { 1.0f, 1.0f });

		CHECK_THROWS(g.set(3, 0, 0, { 0.0f, 0.0f }), std::out_of_range);
		CHECK_THROWS(g.set(0, 2, 0, { 0.0f, 0.0f }), std::out_of_range);
		CHECK_THROWS(g.set(0, 0, 3, { 0.0f, 0.0f }), std::out_of_range);
		CHECK_THROWS(g.set(0, 0, -2, { 0.0f, 0.0f }), std::out_of_range);
		CHECK_THROWS(g.set(0, 0, std::numeric_limits<int>::min(), { 0.0f, 0.0f }), std::out_of_range);
		CHECK_THROWS(g.set(0, 0, 1, { 0.0f }), std::invalid_argument);
		CHECK_THROWS(g.set(0, 0, 1, { 0.0f, 0.0f, 0.0f }), std::invalid_argument);
		CHECK_THROWS(g.set_weights(0, 0, {}), std::invalid_argument);
		CHECK_THROWS(g.set_weights(2, 5, { 0.0f, 0.0f }), std::out_of_range);
		CHECK_THROWS(g.at(3, 0), std::out_of_range);
		CHECK_THROWS(g.at(0, 2), std::out_of_range);
		const PackedAdjacency& cg = g;
		CHECK_THROWS(cg.at(0, 2), std::out_of_range);

		// Failed writes leave the slot unchanged.
		CHECK(g(0, 0).target == 1 && g(0, 0).weight[0] == 1.0f && g(0, 0).weight[1] == 1.0f);

		// Bad dimensions: graph unchanged afterwards.
		const std::size_t max_size = std::numeric_limits<std::size_t>::max();
		const std::size_t too_many_rows = static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()) + 1;
		CHECK_THROWS(g.reset(too_many_rows, 1, 1), std::length_error);
		CHECK_THROWS(g.reset(1u << 20, max_size / 2, 1), std::length_error);
		CHECK_THROWS(g.reset(1u << 16, 1u << 16, max_size / 2), std::length_error);
		CHECK_THROWS(g.reset(1, max_size / 2, 1), std::length_error); // exceeds vector max_size
		CHECK(g.size() == 3 && g.cols() == 2 && g.channels() == 2);
		CHECK(g(0, 0).target == 1 && g(0, 0).weight[1] == 1.0f);
	}

	void test_storage_bytes() {
		PackedAdjacency g;
		g.reset(100, 12, 3);
		const std::size_t expected = 100 * 12 * (sizeof(std::int32_t) + 3 * sizeof(float));
		CHECK(g.storage_bytes() == expected);
		CHECK(g.storage_bytes() == 100 * 12 * 16);
	}

	void test_copy_move() {
		PackedAdjacency a;
		a.reset(3, 2, 2);
		a.set(2, 1, 0, { 4.0f, 8.0f });

		PackedAdjacency b(a);
		CHECK(b.size() == 3 && b.cols() == 2 && b.channels() == 2);
		CHECK(b(2, 1).target == 0 && b(2, 1).weight[1] == 8.0f);
		b(2, 1).weight[1] = 1.0f;
		b(2, 1).target = 1;
		CHECK(a(2, 1).weight[1] == 8.0f && a(2, 1).target == 0); // deep copy
		CHECK(b(2, 1).weight.data() != a(2, 1).weight.data());

		PackedAdjacency c;
		c.reset(1, 1, 1);
		c = a;
		CHECK(c.size() == 3 && c(2, 1).weight[0] == 4.0f);

		const float* a_ptr = a(2, 1).weight.data();
		PackedAdjacency d(std::move(a));
		CHECK(d.size() == 3 && d(2, 1).weight[1] == 8.0f);
		CHECK(d(2, 1).weight.data() == a_ptr); // buffer ownership transferred
		CHECK(a.size() == 0 && a.cols() == 0 && a.channels() == 0 && a.storage_bytes() == 0);

		PackedAdjacency e;
		e.reset(2, 2, 2);
		e = std::move(d);
		CHECK(e.size() == 3 && e(2, 1).target == 0);
		CHECK(d.size() == 0 && d.storage_bytes() == 0);

		// Moved-from object is reusable.
		d.reset(1, 1, 1);
		d.set(0, 0, 0, { 3.0f });
		CHECK(d(0, 0).weight[0] == 3.0f);

		// Self-assignment.
		PackedAdjacency& e_ref = e;
		e = e_ref;
		CHECK(e.size() == 3 && e(2, 1).weight[0] == 4.0f);
		e = std::move(e_ref);
		CHECK(e.size() == 3 && e(2, 1).weight[0] == 4.0f);
	}

} // namespace

int main() {
	test_empty();
	test_multichannel_read_write();
	test_zero_channels();
	test_reset_shape();
	test_bad_inputs();
	test_storage_bytes();
	test_copy_move();

	if (g_failures != 0) {
		std::fprintf(stderr, "packed_adjacency_test: %d failure(s)\n", g_failures);
		return 1;
	}
	std::printf("packed_adjacency_test: all checks passed\n");
	return 0;
}
