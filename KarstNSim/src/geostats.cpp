/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr for new methods + modifications to original methods
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

Rewritten entirely (except the variogram_value fonction, written by G. Rongier in 2015) from the pseudo-code algorithm (Algorithm 2) in the paper of Frantz et al., 2021 "Analysis and stochastic simulation of geometrical properties of conduits in karstic networks".
This work was performed in the frame of the RING project at Université de Lorraine.

***************************************************************/

#include "KarstNSim/geostats.h"

std::vector<int> find_neighborhood(
	int current_node_index,
	const KarstNSim::KarsticSkeleton* const curve,
	int number_max_of_neighborhood_points,
	float range_of_neighborhood,
	std::string type_neighborhood,
	const std::vector<float>& node_values
) {
	// --- Preconditions & early exits ------------------------------------------------------------
	const int N = int(curve->nodes.size());
	std::vector<int> neighborhood;

	if (N == 0) return neighborhood;
	if (current_node_index < 0 || current_node_index >= N) return neighborhood;
	if (number_max_of_neighborhood_points <= 0 || range_of_neighborhood <= 0.f) return neighborhood;

	const bool has_node_values = (node_values.size() == size_t(N));

	// Branch constraint (for "branch" neighborhood type)
	const int branch_id = curve->nodes.at(current_node_index).branch_id;
	const bool restrict_to_branch = (type_neighborhood == "branch");

	// --- Dijkstra truncated by range & capped by K ----------------------------------------------
	struct NodeKey {
		float dist;
		int idx;
		bool operator>(const NodeKey& other) const { return dist > other.dist; }
	};

	const float INF = std::numeric_limits<float>::infinity();
	std::vector<float> dist(N, INF);
	std::vector<char>  visited(N, 0);

	auto edge_length = [&](int u, int v) -> float {
		// Geometric edge weight; same logic as before (Euclidean length).
		const Vector3& pu = curve->nodes[u].p;
		const Vector3& pv = curve->nodes[v].p;
		return KarstNSim::magnitude(pu - pv);
	};

	std::priority_queue<NodeKey, std::vector<NodeKey>, std::greater<NodeKey>> pq;
	dist[current_node_index] = 0.f;
	pq.push({ 0.f, current_node_index });

	// We will collect candidates in (idx, dist) to sort and truncate deterministically.
	std::vector<std::pair<int, float>> candidates;
	candidates.reserve(std::min<int>(number_max_of_neighborhood_points * 2, 128));

	// Safety cap: avoids pathological blow-ups in adversarial graphs
	const size_t MAX_EXPANSIONS = std::max<size_t>(1000000, size_t(N) * 8);
	size_t popped = 0;

	while (!pq.empty()) {
		const auto top = pq.top();
		pq.pop();
		const float du = top.dist;
		const int u = top.idx;
		++popped;

		if (popped > MAX_EXPANSIONS) break;
		if (u < 0 || u >= N) continue;
		if (visited[u]) continue;
		visited[u] = 1;

		// Range truncation: in Dijkstra, if the minimum key exceeds R, remaining keys are ≥ du > R.
		if (du > range_of_neighborhood) break;

		// Skip the seed itself; accept valid neighbors meeting filters.
		if (u != current_node_index) {
			bool accept = true;

			if (restrict_to_branch && curve->nodes[u].branch_id != branch_id)
				accept = false;

			if (accept && has_node_values) {
				// NDV filter used previously: node_values[i] == -99999 means "no data"
				const float val = node_values[u];
				const bool not_ndv = std::fabs(val - (-99999.f)) > 1e-12f;
				if (!not_ndv) accept = false;
			}

			if (accept)
				candidates.emplace_back(u, du);
		}

		// Explore neighbors
		const auto& conns = curve->nodes[u].connections;
		for (const auto& conn : conns) {
			const int v = conn.destindex;
			if (v < 0 || v >= N) continue;
			if (visited[v]) continue;

			const float w = edge_length(u, v);
			if (!std::isfinite(w) || w < 0.f) continue;

			const float alt = du + w;
			if (!std::isfinite(alt) || alt > range_of_neighborhood) continue;

			if (alt < dist[v]) {
				dist[v] = alt;
				pq.push({ alt, v });
			}
		}
	}

	// Sort by metric distance and truncate to K nearest within range
	std::sort(candidates.begin(), candidates.end(),
		[](const auto& a, const auto& b) { return a.second < b.second; });

	if (int(candidates.size()) > number_max_of_neighborhood_points)
		candidates.resize(number_max_of_neighborhood_points);

	neighborhood.reserve(candidates.size());
	for (const auto& kv : candidates)
		neighborhood.push_back(kv.first);

	return neighborhood;
}

// variogram_value (Rongier-2015)
float variogram_value(
	const float& distance,
	const float& sill,
	const float& nugget,
	const float& range,
	const std::string& vario_model) {

	float vario_value = -1;
	// Check the input parameters.

	if (nugget >= 0. && sill >= nugget && range >= 0. && vario_model != "") {
		// Get the variogram value following the model of variogram.
		if (vario_model == "Gaussian") {
			vario_value =
				nugget +
				(sill - nugget) *
				(1 - exp(-(3 * distance / range) * (distance / range)));
		}
		else if (vario_model == "Spherical") {
			if (distance <= range) {
				vario_value =
					nugget +
					(sill - nugget) *
					((3 / 2) * (distance / range) - (1 / 2) * pow((distance / range), 3.));
			}
			else if (distance > range) {
				vario_value = sill;
			}
		}
		else if (vario_model == "Exponential") {
			vario_value =
				nugget +
				(sill - nugget) *
				(1 - exp(-3 * distance / range));
		}
		else if (vario_model == "Nugget") {
			vario_value = nugget;
		}
	}
	return vario_value;
}

// Function to calculate the inverse of a square matrix using Gaussian elimination
bool invert_matrix(const std::vector<std::vector<float>>& input, std::vector<std::vector<float>>& output) {
	// Check if the input matrix is square
	int size = int(input.size());
	if (size != input[0].size()) {
		std::cerr << "Input matrix is not square." << std::endl;
		return false;
	}

	// Initialize the output matrix as the identity matrix
	output = std::vector<std::vector<float>>(size, std::vector<float>(size, 0.0));
	for (int i = 0; i < size; ++i) {
		output[i][i] = 1.0;
	}

	// Copy the input matrix to avoid modifying the original
	std::vector<std::vector<float>> A = input;

	// Gaussian elimination with partial pivoting
	for (int i = 0; i < size; ++i) {
		// Find the pivot row
		int max_row = i;
		for (int k = i + 1; k < size; ++k) {
			if (std::abs(A[k][i]) > std::abs(A[max_row][i])) {
				max_row = k;
			}
		}

		// Swap the current row with the pivot row
		if (max_row != i) {
			A[i].swap(A[max_row]);
			output[i].swap(output[max_row]);
		}

		// Check if the matrix is singular
		if (A[i][i] == 0.0) {
			std::cerr << "Matrix is singular." << std::endl;
			return false;
		}

		// Scale the current row to make the diagonal element 1
		float scale = 1.0 / A[i][i];
		for (int j = 0; j < size; ++j) {
			A[i][j] *= scale;
			output[i][j] *= scale;
		}

		// Eliminate the current column
		for (int k = 0; k < size; ++k) {
			if (k != i) {
				float factor = A[k][i];
				for (int j = 0; j < size; ++j) {
					A[k][j] -= factor * A[i][j];
					output[k][j] -= factor * output[i][j];
				}
			}
		}
	}
	return true;
}

// --- Helper: truncated Dijkstra shortest-path distances to a set of targets ---
// Returns distances from 'src' to each 'targets[t]' (INF if beyond 'range_cap' or unreachable).
std::vector<float> dijkstra_to_targets_truncated(
	int src,
	const KarstNSim::KarsticSkeleton* curve,
	const std::vector<int>& targets,
	float range_cap
) {
	const int N = (int)curve->nodes.size();
	const float INF = std::numeric_limits<float>::infinity();

	std::vector<float> dist(N, INF);
	std::vector<char>  seen(N, 0);

	struct Q { float d; int i; bool operator>(const Q& o) const { return d > o.d; } };
	std::priority_queue<Q, std::vector<Q>, std::greater<Q>> pq;

	auto edge_len = [&](int u, int v)->float {
		const Vector3& pu = curve->nodes[u].p;
		const Vector3& pv = curve->nodes[v].p;
		return KarstNSim::magnitude(pu - pv);
	};

	dist[src] = 0.f;
	pq.push({ 0.f, src });

	// For an early-out when every target is settled
	std::vector<char> is_target(N, 0);
	for (int t : targets) if (t >= 0 && t < N) is_target[t] = 1;
	int remaining = (int)targets.size();

	while (!pq.empty()) {
		auto [du, u] = pq.top(); pq.pop();
		if (seen[u]) continue;
		seen[u] = 1;

		if (du > range_cap) break;               // all further keys ≥ du
		if (is_target[u]) {                      // settled a target
			if (--remaining == 0) break;        // all found
		}

		for (const auto& c : curve->nodes[u].connections) {
			const int v = c.destindex;
			if (v < 0 || v >= N) continue;
			if (seen[v]) continue;
			const float w = edge_len(u, v);
			if (!std::isfinite(w) || w < 0.f) continue;
			const float alt = du + w;
			if (alt >= dist[v] || alt > range_cap) continue;
			dist[v] = alt;
			pq.push({ alt, v });
		}
	}

	std::vector<float> out; out.reserve(targets.size());
	for (int t : targets) out.push_back((t >= 0 && t < N) ? dist[t] : INF);
	return out;
}

// --- Kriging that builds a local distance matrix on the fly (no global distance_mat) ---
void kriging_in_point_on_the_fly(
	int current_node_index,
	const std::vector<int>& neighborhood,          // indices of conditioning nodes
	const std::vector<float>* kriging_distribution,
	const KarstNSim::KarsticSkeleton* curve,
	const float& vario_range,
	const float& vario_sill,
	const float& vario_nugget,
	const std::string& vario_model,
	const std::vector<float>& node_values,         // gaussian values incl. -99999 NDV
	float& var_estimation,
	float& val_estimation,
	float range_cap                                  // SAME radius used to build the neighborhood
) {
	// If the neighborhood is empty, sample directly (preserves previous behavior).
	if (neighborhood.empty()) {
		val_estimation = select_random_element(*kriging_distribution);
		// Provide a nominal variance (sill) so downstream logs don’t see garbage.
		var_estimation = vario_sill;
		return;
	}

	const int K = (int)neighborhood.size();
	// 1) Distances current -> neighbors (shortest path, truncated at range_cap)
	std::vector<int> targets = neighborhood;
	std::vector<float> d_cur_to_nb = dijkstra_to_targets_truncated(
		current_node_index, curve, targets, range_cap
	);

	// 2) Pairwise neighbor-neighbor distances: run a truncated Dijkstra from each neighbor,
	//    but only extract distances to the other neighbors (K is small: ≤ ~16).
	std::vector<std::vector<float>> d_nb_to_nb(K, std::vector<float>(K, 0.f));
	for (int i = 0; i < K; ++i) {
		std::vector<float> di = dijkstra_to_targets_truncated(neighborhood[i], curve, targets, range_cap);
		for (int j = 0; j < K; ++j) d_nb_to_nb[i][j] = di[j];
	}

	// 3) Build local variogram/covariance system (Kriging with Lagrange multiplier)
	//    Matrix size = (K+1) x (K+1)
	const int M = K + 1;
	std::vector<std::vector<float>> Gamma(M, std::vector<float>(M, 0.f));
	std::vector<float>            rhs(M, 0.f);

	auto gamma = [&](float h)->float {
		if (h <= 0.f) return 0.f;
		return variogram_value(h, vario_sill, vario_nugget, vario_range, vario_model);
	};

	// Fill Γ (0..K-1, 0..K-1) with γ(h_ij)
	for (int i = 0; i < K; ++i) {
		for (int j = 0; j < K; ++j) {
			float hij = d_nb_to_nb[i][j];
			// If unreachable within range_cap, set γ ≈ sill (i.e., no correlation).
			Gamma[i][j] = std::isfinite(hij) ? gamma(hij) : vario_sill;
		}
	}
	// Lagrange row/col
	for (int i = 0; i < K; ++i) {
		Gamma[i][K] = 1.f;
		Gamma[K][i] = 1.f;
	}
	Gamma[K][K] = 0.f;

	// Right-hand side: γ(h_i0) where 0 is "current node"
	for (int i = 0; i < K; ++i) {
		float hi0 = d_cur_to_nb[i];
		rhs[i] = std::isfinite(hi0) ? gamma(hi0) : vario_sill;
	}
	rhs[K] = 1.f;

	// Solve Γ * [λ; μ] = rhs
	std::vector<std::vector<float>> Gamma_inv;
	if (!invert_matrix(Gamma, Gamma_inv)) {
		// Numerical fallback: if inversion fails, sample from distribution
		var_estimation = vario_sill;
		val_estimation = select_random_element(*kriging_distribution);
		return;
	}

	std::vector<float> w(M, 0.f);
	for (int i = 0; i < M; ++i) {
		float acc = 0.f;
		for (int j = 0; j < M; ++j) acc += Gamma_inv[i][j] * rhs[j];
		w[i] = acc;
	}

	// Kriging estimate at 0: z* = sum_i λ_i * z_i
	float estimate = 0.f;
	for (int i = 0; i < K; ++i) estimate += w[i] * node_values[neighborhood[i]];
	val_estimation = estimate;

	// Kriging variance: σ² = γ(0) - sum_i λ_i * γ(h_i0) - μ  (γ(0) = 0 by definition)
	float cross = 0.f;
	for (int i = 0; i < K; ++i) cross += w[i] * rhs[i];
	float mu = w[K];
	var_estimation = std::max(0.f, -cross - mu);
}

void kriging_in_point(
	const int& current_node_index,
	const std::vector<int>& neighborhood,
	const std::vector<float>* kriging_distribution,
	const Array2D<float>& mat_distance,
	const float& vario_range,
	const float& vario_sill,
	const float& vario_nugget,
	const std::string& vario_model,
	const std::vector<float>& node_values, // node values, not updated here
	float& var_estimation, // variance estimated, used for SGS simulation
	float& val_estimation // value estimated, used for SGS simulation
) {

	float average = std::accumulate(kriging_distribution->begin(), kriging_distribution->end(), 0.0) / kriging_distribution->size();

	val_estimation = -99999;

	// If the neighborhood is empty, perform a random draw in the initial distribution.
	if (neighborhood.empty()) {
		val_estimation = select_random_element(*kriging_distribution);
		return;
	}
	else {
		int nb_neigh = int(neighborhood.size());

		// Create matrix of distances that can be used for calculation
		std::vector<std::vector<float>> local_mat_distance(nb_neigh + 1, std::vector<float>(nb_neigh + 1));
		std::vector<std::vector<float>> local_mat_vario_values(nb_neigh + 1, std::vector<float>(nb_neigh + 1));
		std::vector<std::vector<float>> local_mat_cov(nb_neigh + 1, std::vector<float>(nb_neigh + 1));

		local_mat_distance[0][0] = 0.;
		for (int i = 0; i < nb_neigh; ++i) {
			local_mat_distance[i + 1][0] = mat_distance(current_node_index, neighborhood[i]);
			local_mat_distance[0][i + 1] = mat_distance(current_node_index, neighborhood[i]);
		}

		// Fill the rest of the matrix
		for (int i = 0; i < nb_neigh; ++i) {
			for (int j = 0; j < nb_neigh; ++j) {
				local_mat_distance[i + 1][j + 1] = mat_distance(neighborhood[i], neighborhood[j]);
			}
		}

		float a = vario_range;
		float C = vario_sill;
		float C0 = vario_nugget;

		// Fill matrix of variogram values
		for (int j = 0; j < nb_neigh + 1; ++j) {
			float h = local_mat_distance[0][j];
			if (h == 0.) { // It should not be the case
				local_mat_vario_values[0][j] = 0.;
			}
			else {
				local_mat_vario_values[0][j] = variogram_value(h, C, C0, a, vario_model);
			}
			local_mat_vario_values[j][0] = local_mat_vario_values[0][j];
		}

		for (int i = 1; i < nb_neigh + 1; ++i) {
			for (int j = i; j < nb_neigh + 1; ++j) {
				float h = local_mat_distance[i][j]; // Matrix of values in variogram
				if (h == 0.) {
					local_mat_vario_values[i][j] = 0.;
				}
				else {
					local_mat_vario_values[i][j] = variogram_value(h, C, C0, a, vario_model);
				}
				local_mat_vario_values[j][i] = local_mat_vario_values[i][j];
			}
		}

		// Matrix of covariance
		for (int i = 0; i < nb_neigh + 1; ++i) {
			for (int j = 0; j < nb_neigh + 1; ++j) {
				local_mat_cov[i][j] = C - local_mat_vario_values[i][j];
			}
		}

		// Perform simple kriging
		std::vector<std::vector<float>> k0(nb_neigh, std::vector<float>(1));

		for (int i = 0; i < nb_neigh; ++i) {
			k0[i][0] = local_mat_cov[0][i + 1];
		}

		std::vector<std::vector<float>> K(nb_neigh, std::vector<float>(nb_neigh));

		for (int i = 0; i < nb_neigh; ++i) {
			for (int j = 0; j < nb_neigh; ++j) {
				K[i][j] = local_mat_cov[i + 1][j + 1];
			}
		}

		std::vector<std::vector<float>> lambda(nb_neigh, std::vector<float>(1));
		std::vector<std::vector<float>> inv_K(nb_neigh, std::vector<float>(nb_neigh));

		bool invert_success = invert_matrix(K, inv_K);

		if (!invert_success) {
			// Handle inversion failure
				return;
		}

		for (int i = 0; i < nb_neigh; ++i) {
			lambda[i][0] = 0;
			for (int j = 0; j < nb_neigh; ++j) {
				lambda[i][0] += inv_K[i][j] * k0[j][0];
			}
		}

		val_estimation = average;
		float sum_lambda = 0;

		for (int i = 0; i < nb_neigh; ++i) {
			val_estimation += lambda[i][0] * node_values[neighborhood[i]];
			sum_lambda += lambda[i][0];
		}

		val_estimation -= sum_lambda * average;

		float var_temp = 0;

		for (int i = 0; i < nb_neigh; ++i) {
			var_temp += k0[i][0] * lambda[i][0];
		}

		var_estimation = local_mat_cov[0][0] - var_temp;

		return;
	}
}



void save_data(const std::vector<float>& data, const std::string& filename) {
	std::ofstream file(filename);

	for (const auto& value : data) {
		file << value << std::endl;
	}

	file.close();

}

// ===== Upstream orientation from springs (component-wise) ====================
// The graph is split into connected components; for each component we:
//  1) Map input springs (3D positions) to their nearest skeleton node;
//  2) Keep only the springs that fall inside the component;
//  3) Among those, select the lowest-Z spring and any other spring in the
//     component whose Z <= lowest_Z + 40 m (multi-source rule);
//  4) Run a multi-source Dijkstra, direct edges by strictly increasing distance,
//     then accumulate total upstream edge-union length (each edge counted once).

namespace {
	constexpr float DCURV_EPS = 1e-6f;

	inline float edge_length(const KarstNSim::KarsticSkeleton* sk, int a, int b) {
		return KarstNSim::magnitude(sk->nodes[a].p - sk->nodes[b].p);
	}

	// --- Utilities ------------------------------------------------------------

	inline float sqr(float x) { return x * x; }

	// Undirected edge key (works while N < 2^32).
	inline uint64_t make_edge_key(int a, int b) {
		if (a > b) std::swap(a, b);
		return ((uint64_t)(uint32_t)a << 32) | (uint64_t)(uint32_t)b;
	}

	static std::vector<int> build_components(const KarstNSim::KarsticSkeleton* sk) {
		const int N = (int)sk->nodes.size();
		std::vector<int> comp(N, -1);
		int cid = 0;
		std::vector<int> stack; stack.reserve(N);

		for (int s = 0; s < N; ++s) {
			if (comp[s] != -1) continue;
			comp[s] = cid;
			stack.clear();
			stack.push_back(s);
			while (!stack.empty()) {
				int v = stack.back(); stack.pop_back();
				for (const auto& c : sk->nodes[v].connections) {
					int u = c.destindex;
					if (u >= 0 && u < N && comp[u] == -1) {
						comp[u] = cid;
						stack.push_back(u);
					}
				}
			}
			++cid;
		}
		return comp;
	}

	static int nearest_node_index(const KarstNSim::KarsticSkeleton* sk,
		const Vector3& P)
	{
		// Linear scan is fine here; can be replaced by a spatial index if needed.
		int best = -1;
		float best2 = std::numeric_limits<float>::infinity();
		for (int i = 0; i < (int)sk->nodes.size(); ++i) {
			const auto& Q = sk->nodes[i].p;
			float d2 = sqr(P.x - Q.x) + sqr(P.y - Q.y) + sqr(P.z - Q.z);
			if (d2 < best2) { best2 = d2; best = i; }
		}
		return best;
	}

	// Select sources per component using the *input* spring Z, then map to nodes.
	// We first keep all input springs whose nearest node lies in the component,
	// compute zmin on the *spring Z*, and keep only those with S.z <= zmin + z_window.
	// Finally we map these selected springs to their nearest node indices and deduplicate.
	static std::vector<int> pick_component_sources_from_springs(
		const KarstNSim::KarsticSkeleton* sk,
		const std::vector<int>& comp,
		int cid,
		const std::vector<Vector3>& springs_xyz,
		float z_window = 40.f)
	{
		struct Cand { int spring_id; int node_id; float spring_z; };
		std::vector<Cand> cand;
		cand.reserve(springs_xyz.size());

		// Map every spring to its nearest node and keep only those landing in this component
		for (int sid = 0; sid < (int)springs_xyz.size(); ++sid) {
			const auto& S = springs_xyz[sid];
			const int nidx = nearest_node_index(sk, S);
			if (nidx >= 0 && comp[nidx] == cid) {
				cand.push_back({ sid, nidx, S.z }); // NOTE: keep *spring* Z here
			}
		}
		if (cand.empty()) return {};

		// Lowest *spring* Z among candidates of this component
		float zmin = cand[0].spring_z;
		for (const auto& c : cand) zmin = std::min(zmin, c.spring_z);

		// Keep only springs within +z_window of the lowest spring Z (component-wise)
		std::vector<int> sources_nodes;
		sources_nodes.reserve(cand.size());
		for (const auto& c : cand) {
			if (c.spring_z <= zmin + z_window) {
				sources_nodes.push_back(c.node_id);
			}
		}

		// Deduplicate node indices (several springs could map to the same node)
		std::sort(sources_nodes.begin(), sources_nodes.end());
		sources_nodes.erase(std::unique(sources_nodes.begin(), sources_nodes.end()), sources_nodes.end());

		return sources_nodes;
	}

	// --- Multi-source shortest paths (Dijkstra) with source label propagation ----

	// Overload that can optionally return the "nearest source label" for each node.
	// If label_out != nullptr, label_out->size() will be set to N with:
	//   -1 = unreachable/different component; otherwise index in `sources` (0..S-1).
	static std::vector<float> dijkstra_multi_source(
		const KarstNSim::KarsticSkeleton* sk,
		const std::vector<int>& comp, int target_cid,
		const std::vector<int>& sources,
		std::vector<int>* label_out)
	{
		const int N = (int)sk->nodes.size();
		const float INF = std::numeric_limits<float>::infinity();
		std::vector<float> dist(N, INF);
		if (label_out) label_out->assign(N, -1);

		struct Item { float d; int v; int src_id; };
		struct Cmp { bool operator()(const Item& a, const Item& b) const { return a.d > b.d; } };
		std::priority_queue<Item, std::vector<Item>, Cmp> pq;

		// Initialize each selected spring as an independent source with its own label
		for (int si = 0; si < (int)sources.size(); ++si) {
			const int s = sources[si];
			dist[s] = 0.0f;
			if (label_out) (*label_out)[s] = si;
			pq.push({ 0.0f, s, si });
		}

		while (!pq.empty()) {
			auto cur = pq.top(); pq.pop();
			const float dv = cur.d; const int v = cur.v; const int vlabel = cur.src_id;
			if (dv > dist[v]) continue;
			if (comp[v] != target_cid) continue;

			for (const auto& c : sk->nodes[v].connections) {
				const int u = c.destindex;
				if (u < 0 || u >= N) continue;
				if (comp[u] != target_cid) continue;

				const float w = KarstNSim::magnitude(sk->nodes[v].p - sk->nodes[u].p);
				const float nd = dv + w;

				if (nd + DCURV_EPS < dist[u]) {
					dist[u] = nd;
					if (label_out) (*label_out)[u] = vlabel;
					pq.push({ nd, u, vlabel });
				}
				// Optional: tie-break on equal distances -> smallest label for determinism
				else if (label_out && std::abs(nd - dist[u]) <= DCURV_EPS) {
					if (vlabel < (*label_out)[u]) {
						(*label_out)[u] = vlabel;
						// No need to push again: distance didn't change.
					}
				}
			}
		}
		return dist;
	}

	/// @brief Compute per-basin edge-length totals using a multi-source Voronoi
///        partition on the graph. For each undirected edge (u,v):
///        - If labels match: add full length to that basin;
///        - If labels differ: split the edge at the equidistance point:
///              du + x = dv + (w - x) -> x = (dv - du + w)/2,
///          and add portions x and (w-x) to the two basins respectively.
///        This guarantees that the sum over basins equals the component's total
///        undirected edge length.
/// @param curve     Skeleton
/// @param comp,cid  Component labelling & current component id
/// @param dist      Dijkstra distances (0 at sources, increasing upstream)
/// @param label_of  For each node, index in `sources` of nearest source (0..S-1), or -1
/// @param sources   List of source node indices (as passed to Dijkstra)
/// @param basin_len_out Filled with S entries: total length per basin
/// @param total_len_out  Filled with component total undirected edge length
	static void accumulate_basin_edge_lengths(
		const KarstNSim::KarsticSkeleton* curve,
		const std::vector<int>& comp, int cid,
		const std::vector<float>& dist,
		const std::vector<int>& label_of,
		const std::vector<int>& sources,
		std::vector<double>& basin_len_out,
		double& total_len_out)
	{
		const int N = (int)curve->nodes.size();
		const int S = (int)sources.size();
		basin_len_out.assign(S, 0.0);
		total_len_out = 0.0;

		for (int v = 0; v < N; ++v) if (comp[v] == cid) {
			for (const auto& c : curve->nodes[v].connections) {
				const int u = c.destindex;
				if (u < 0 || u >= N || comp[u] != cid) continue;
				if (u >= v) continue; // undirected: count each edge once

				const float w = KarstNSim::magnitude(curve->nodes[u].p - curve->nodes[v].p);
				total_len_out += w;

				const int lu = label_of[u];
				const int lv = label_of[v];

				// If either endpoint has no label (should not happen inside cid), skip split
				if (lu < 0 || lv < 0) continue;

				if (lu == lv) {
					basin_len_out[lu] += w;
				}
				else {
					// Split the edge proportionally at equidistance along the segment.
					// Solve du + x = dv + (w - x) => x = (dv - du + w) / 2
					const float du = dist[u];
					const float dv = dist[v];
					float x = (dv - du + w) * 0.5f;
					if (x < 0.0f)       x = 0.0f;
					if (x > w)          x = w;
					basin_len_out[lu] += (double)x;
					basin_len_out[lv] += (double)(w - x);
				}
			}
		}
	}

	/// @brief Save the "basin boundary" information to CSV for 2D plotting.
	///			To be used for debugging purposes
	///        Each row = one graph edge (u,v). If labels differ, we emit the
	///        equidistance point (split_x, split_y); otherwise fields are empty.
	///        Columns:
	///        u,v, lu,lv,  ux,uy, vx,vy,  split_x,split_y,  w
	static void save_basin_edges_csv(
		const KarstNSim::KarsticSkeleton* curve,
		const std::vector<int>& comp, int cid,
		const std::vector<float>& dist,
		const std::vector<int>& label_of,
		const std::string& out_csv_path)
	{
		std::ofstream ofs(out_csv_path);
		if (!ofs) {
			return;
		}
		ofs << "u,v,lu,lv,ux,uy,vx,vy,split_x,split_y,w\n";

		const int N = (int)curve->nodes.size();
		for (int v = 0; v < N; ++v) if (comp[v] == cid) {
			for (const auto& c : curve->nodes[v].connections) {
				const int u = c.destindex;
				if (u < 0 || u >= N || comp[u] != cid) continue;
				if (u >= v) continue; // undirected: once

				const auto& Pu = curve->nodes[u].p;
				const auto& Pv = curve->nodes[v].p;
				const float w = KarstNSim::magnitude(Pu - Pv);
				const int lu = label_of[u];
				const int lv = label_of[v];

				if (lu != -1 && lv != -1 && lu != lv) {
					// Equidistance point along segment (2D plan view: x,y)
					const float du = dist[u];
					const float dv = dist[v];
					float x = (dv - du + w) * 0.5f;
					if (x < 0.0f)       x = 0.0f;
					if (x > w)          x = w;
					const float t = (w > 0.0f) ? (x / w) : 0.5f;
					const float sx = Pu.x + t * (Pv.x - Pu.x);
					const float sy = Pu.y + t * (Pv.y - Pu.y);
					ofs << u << ',' << v << ',' << lu << ',' << lv << ','
						<< Pu.x << ',' << Pu.y << ','
						<< Pv.x << ',' << Pv.y << ','
						<< sx << ',' << sy << ','
						<< w << '\n';
				}
				else {
					ofs << u << ',' << v << ',' << lu << ',' << lv << ','
						<< Pu.x << ',' << Pu.y << ','
						<< Pv.x << ',' << Pv.y << ','
						<< ',' << ',' << ','   // empty split_x,split_y,w if unlabeled? keep w
						<< w << '\n';
				}
			}
		}
		ofs.close();
	}

	// --- Build upstream DAG (multi-parents) ----------------------------------

	static void build_upstream_parents(
		const KarstNSim::KarsticSkeleton* sk,
		const std::vector<int>& comp, int cid,
		const std::vector<float>& dist,
		std::vector<std::vector<int>>& parents)
	{
		const int N = (int)sk->nodes.size();
		parents.assign(N, {});
		for (int v = 0; v < N; ++v) {
			if (comp[v] != cid) continue;
			for (const auto& c : sk->nodes[v].connections) {
				int u = c.destindex;
				if (u < 0 || u >= N) continue;
				if (comp[u] != cid) continue;
				if (dist[u] > dist[v] + DCURV_EPS) {
					parents[v].push_back(u); // u is strictly upstream of v
				}
			}
		}
	}

	// --- Exact upstream union per target set ---------------------------------
	// Each edge is counted at most once per target by using an unordered_set of undirected edge keys.

	static std::vector<float> exact_dcurv_for_targets(
		const KarstNSim::KarsticSkeleton* curve,
		const std::vector<int>& comp, int cid,
		const std::vector<std::vector<int>>& parents,
		const std::vector<int>& targets)
	{
		using EdgeKey = uint64_t;
		std::vector<float> out; out.reserve(targets.size());

		std::vector<int> stack; stack.reserve(1024);

		for (int t : targets) {
			if (t < 0 || t >= (int)curve->nodes.size() || comp[t] != cid) {
				out.push_back(0.0f);
				continue;
			}

			std::unordered_set<EdgeKey> seen_edges;
			seen_edges.reserve(1024); // heuristic

			float total = 0.0f;
			stack.clear();
			stack.push_back(t);

			// Non-recursive DFS up the parents DAG
			while (!stack.empty()) {
				const int v = stack.back();
				stack.pop_back();

				for (int u : parents[v]) {
					const EdgeKey ek = make_edge_key(u, v);
					if (seen_edges.insert(ek).second) {
						total += edge_length(curve, u, v);
						stack.push_back(u);
					}
					// else: edge already counted; skip to avoid double counting.
				}
			}

			out.push_back(total);
		}

		return out;
	}
} // anonymous namespace


/// @brief Compute total upstream curvilinear length per node, using only springs
///        provided by the caller and handling disconnected components independently.
/// @param curve         Skeleton graph
/// @param springs_xyz   3D positions of springs (as provided by KarsticNetwork::pt_spring)
/// @return              Vector dcurv: for each node v, sum of lengths of all edges
///                      in the upstream subgraph that feeds v (component-wise, edge-union).
std::vector<float> compute_upstream_curvilinear_length(
	const KarstNSim::KarsticSkeleton* curve,
	const std::vector<Vector3>& springs_xyz)
{
	const int N = (int)curve->nodes.size();
	std::vector<float> dcurv(N, 0.0f);

	// 1) Connected components
	std::vector<int> comp = build_components(curve);
	int ncomp = 0; for (int x : comp) ncomp = std::max(ncomp, x + 1);

	// 2) Process each component independently
	for (int cid = 0; cid < ncomp; ++cid) {
		// 2.1) Select sources among springs located in this component:
		//      lowest-Z spring + any other spring within +40 m in Z.
		std::vector<int> sources = pick_component_sources_from_springs(
			curve, comp, cid, springs_xyz, 40.f);

		if (sources.empty()) {
			//	"dcurv will be 0 for all nodes in this component.");
			continue; // nothing to orient here
		}

		// 2.2) Multi-source Dijkstra distances
		std::vector<float> dist = dijkstra_multi_source(curve, comp, cid, sources, nullptr);

		// 2.2) Multi-source Dijkstra distances + labels (which spring each node belongs to)
		//std::vector<int> label_of; // size N after call, values in [0..S-1] or -1
		//std::vector<float> dist = dijkstra_multi_source(curve, comp, cid, sources, &label_of);

		// --- Basin audit & CSV export ----------------------------------------------
		//std::vector<double> basin_len;
		//double comp_len_check = 0.0;
		//accumulate_basin_edge_lengths(curve, comp, cid, dist, label_of, sources, basin_len, comp_len_check);
		//// Save CSV for plotting boundaries (plan view)
		//const std::string out_csv = "C:\\Users\\gouy2\\PycharmProjects\\project_karst\\outputs\\basin_boundaries_comp" + std::to_string(cid) + ".csv";
		//save_basin_edges_csv(curve, comp, cid, dist, label_of, out_csv);

		// 2.3) Upstream DAG (multi-parents)
		std::vector<std::vector<int>> parents;
		build_upstream_parents(curve, comp, cid, dist, parents);

		// (Optional) Compute total undirected edge length of this component for sanity logs
		double comp_total_len = 0.0;
		for (int v = 0; v < N; ++v) if (comp[v] == cid) {
			for (const auto& c : curve->nodes[v].connections) {
				int u = c.destindex;
				if (u < 0 || u >= N || comp[u] != cid) continue;
				if (u < v) comp_total_len += edge_length(curve, u, v);
			}
		}

		// 2.4) Exact upstream edge-union length for ALL nodes in the component
		std::vector<int> targets; targets.reserve(1024);
		for (int v = 0; v < N; ++v) if (comp[v] == cid) targets.push_back(v);

		std::vector<float> dcurv_comp = exact_dcurv_for_targets(
			curve, comp, cid, parents, targets);

		// Scatter back to global dcurv
		size_t k = 0;
		for (int v = 0; v < N; ++v) {
			if (comp[v] != cid) continue;
			dcurv[v] = dcurv_comp[k++];
			// Optional clamp & warning (should never trigger with union-logic)
			if (dcurv[v] > comp_total_len + 1e-3f) {
				dcurv[v] = (float)comp_total_len;
			}
		}

	}

	// 3) Stats/log
	float mn = std::numeric_limits<float>::infinity();
	float mx = -std::numeric_limits<float>::infinity();
	int n_zero = 0, n_iso = 0;
	for (int i = 0; i < N; ++i) {
		mn = std::min(mn, dcurv[i]);
		mx = std::max(mx, dcurv[i]);
		if (dcurv[i] <= 1e-6f) ++n_zero;
		if (curve->nodes[i].connections.empty()) ++n_iso;
	}

	return dcurv;
}


// === External drift by redundancy-aware weighted OLS with hard-trim outlier rejection ===
// - Predictors:
//     * zwt: vertical distance above the phreatic reference (clamped to 0 below WT)
//     * dcurv: total upstream curvilinear length
//   Both are normalized to [0,1] on the fitting subset.
// - Redundancy weights:
//     * 1D binning along the main active axis (priority: dcurv, else zwt)
//     * individual weight w_i ∝ 1 / (#points in the bin of i), then renormalized to mean=1
// - Class-balance weighting (NEW):
//     * Within the current fitting subset, split total weight mass 50/50
//       between samples located at spring positions (SPRINGS) and all others (OTHERS).
//     * Uses exact proximity test with tolerance EPS = 1e-5 on (x,y,z),
//       via KarstNSim::magnitude(P - S) <= EPS.
// - Outlier rejection (hard 0/1):
//     * compute residuals r_i on the current fit
//     * robust scale via MAD (sigma = 1.4826 * median(|r_i|))
//     * reject i if |r_i| / sigma > C_CUTOFF (Tukey-like hard threshold)
//     * refit on survivors ONLY and RECOMPUTE redundancy + class-balance weights on survivors
// - Geological sign checks:
//     * β_zwt < 0 expected (larger radii near water table, zwt decreases)
//     * β_dcurv > 0 expected (increasing downstream)
// - Returns:
//     * drift (size N) for all nodes
//     * weights_out per node: redundancy-balanced weight for observed survivors,
//       0 for rejected outliers and non-observed nodes.
// -----------------------------------------------------------------------------
std::vector<float> compute_external_drift(
	const KarstNSim::KarsticSkeleton* curve,
	const std::vector<Vector3>& springs_xyz,
	const float& z_phreatic,
	const bool& use_drift_zwt,
	const bool& use_drift_curv,
	const std::vector<float>& eq_radius_values,
	std::vector<float>& weights_out)
{

	const int N = static_cast<int>(curve->nodes.size());
	std::vector<float> drift(N, 0.0f);

	// ---- 0) Build predictors on all nodes -------------------------------------
	// zwt: vertical distance above phreatic surface (clamped at 0 below WT)
	std::vector<float> zwt(N, 0.0f);
	if (use_drift_zwt) {
		for (int i = 0; i < N; ++i) {
			const float raw = curve->nodes[i].p.z - z_phreatic;
			zwt[i] = (raw > 0.0f ? raw : 0.0f);
		}
	}

	// dcurv: total upstream curvilinear length (component-wise, multi-source springs)
	std::vector<float> dcurv(N, 0.0f);
	if (use_drift_curv) {
		dcurv = compute_upstream_curvilinear_length(curve, springs_xyz);
	}

	// ---- 1) Collect observed nodes (non-NDV eq_radius) ------------------------

	const float RADIUS_MAX_FOR_REGRESSION = 2.1f;  // observations with radius > max are excluded (typically ghost-rocks ; to be modified according to need)

	std::vector<int> valid_indices;
	valid_indices.reserve(N);
	for (int i = 0; i < N; ++i) {
		// Keep only finite, non-NDV radii AND enforce radius <= threshold
		if (i < (int)eq_radius_values.size()) {
			const float r = eq_radius_values[i];
			const bool has_data = (std::abs(r - (-99999.0f)) > 1e-12f);
			if (has_data && r <= RADIUS_MAX_FOR_REGRESSION) {
				valid_indices.push_back(i);
			}
		}
	}

	const int n_obs = static_cast<int>(valid_indices.size());
	if (n_obs == 0 || (!use_drift_zwt && !use_drift_curv)) {
		weights_out.assign(N, 0.0f);
		return drift;
	}

	// ---- Per-observation listing -------------------------
	// For each observed datum, report zwt & dcurv used by the regression. 
	for (int id : valid_indices) {
		const Vector3& P = curve->nodes[id].p;
		const float z_val = use_drift_zwt ? zwt[id] : 0.0f;
		const float d_val = use_drift_curv ? dcurv[id] : 0.0f;
	}

	// ---- Helpers ---------------------------------------------------------------
	auto compute_minmax_on = [&](const std::vector<float>& v, const std::vector<int>& idxs) {
		float vmin = std::numeric_limits<float>::infinity();
		float vmax = -std::numeric_limits<float>::infinity();
		for (int id : idxs) { vmin = std::min(vmin, v[id]); vmax = std::max(vmax, v[id]); }
		if (!std::isfinite(vmin) || !std::isfinite(vmax) || std::abs(vmax - vmin) < 1e-12f) {
			vmin = 0.0f; vmax = 1.0f;
		}
		return std::make_pair(vmin, vmax);
	};

	auto redundancy_weights = [&](const std::vector<int>& idxs, const std::vector<float>& axis01) {
		const int M = static_cast<int>(idxs.size());
		std::vector<float> w(M, 1.0f);
		if (M <= 1) return w;

		std::vector<std::pair<float, int>> order; order.reserve(M);
		for (int k = 0; k < M; ++k) order.emplace_back(axis01[k], k);
		std::sort(order.begin(), order.end());

		const int K = std::max(1, std::min(8, M / 10));
		for (int t = 0; t < M; ++t) {
			const int k0 = order[t].second;
			int kL = std::max(0, t - K), kR = std::min(M - 1, t + K);
			float acc = 0.0f, norm = 0.0f;
			for (int s = kL; s <= kR; ++s) {
				const float d = std::abs(order[s].first - order[t].first);
				const float ker = std::max(0.0f, 1.0f - d); // triangular kernel
				acc += ker; norm += 1.0f;
			}
			const float dens = (norm > 0.0f ? acc / norm : 1.0f);
			w[k0] = 1.0f / (1.0f + dens);
		}
		float wmin = *std::min_element(w.begin(), w.end());
		float wmax = *std::max_element(w.begin(), w.end());
		for (float& wi : w) {
			if (wmax > wmin) wi = 0.2f + 0.8f * (wi - wmin) / (wmax - wmin);
			else wi = 1.0f;
		}
		return w;
	};

	// ---- Weighted OLS on a subset with redundancy + class-balance --------------
	auto fit_on_subset = [&](const std::vector<int>& idxs,
		bool use_zwt, bool use_dcurv,
		std::vector<float>& out_beta,
		std::vector<float>& out_weights,
		float& zwt_min, float& zwt_max,
		float& dcurv_min, float& dcurv_max) -> bool
	{
		if (use_zwt)   std::tie(zwt_min, zwt_max) = compute_minmax_on(zwt, idxs);
		if (use_dcurv) std::tie(dcurv_min, dcurv_max) = compute_minmax_on(dcurv, idxs);

		// Axis used for redundancy weighting (priority: dcurv, else zwt), all in [0,1]
		std::vector<float> axis01; axis01.reserve(idxs.size());
		if (use_dcurv) {
			const float rng = std::max(1e-12f, dcurv_max - dcurv_min);
			for (int id : idxs) axis01.push_back((dcurv[id] - dcurv_min) / rng);
		}
		else if (use_zwt) {
			const float rng = std::max(1e-12f, zwt_max - zwt_min);
			for (int id : idxs) axis01.push_back((zwt[id] - zwt_min) / rng);
		}
		out_weights = redundancy_weights(idxs, axis01); // base redundancy weights

		// === Class-balance 50/50 (SPRINGS vs OTHERS), exact proximity (EPS=1e-2) ===
		const float SPR_EPS = 1e-2f;
		float Wspr = 0.0f, Woth = 0.0f;
		size_t Nspr = 0, Noth = 0;
		std::vector<uint8_t> is_spring(idxs.size(), 0); // <-- FIXED: idxs.size()
		for (size_t k = 0; k < idxs.size(); ++k) {
			const int id = idxs[k];
			const Vector3& p = curve->nodes[id].p;

			bool sflag = false;
			for (size_t si = 0; si < springs_xyz.size(); ++si) {
				const float dist = KarstNSim::magnitude(p - springs_xyz[si]);
				// Verbose diagnostic: log magnitude for each pair tested
				if (dist <= SPR_EPS) { sflag = true; break; }
			}

			is_spring[k] = sflag ? 1u : 0u;
			if (sflag) { Wspr += out_weights[k]; ++Nspr; }
			else { Woth += out_weights[k]; ++Noth; }
		}
		const float Wtot = Wspr + Woth;
		if (Wspr > 0.0f && Woth > 0.0f) {
			const float target = 0.5f * Wtot;
			const float fs = target / Wspr;   // factor applied to springs
			const float fo = target / Woth;   // factor applied to others
			for (size_t k = 0; k < idxs.size(); ++k) {
				out_weights[k] *= (is_spring[k] ? fs : fo);
			}
		}
		// ==========================================================================

		const int n_var = 1 + (use_zwt ? 1 : 0) + (use_dcurv ? 1 : 0);
		std::vector<std::vector<float>> XtX(n_var, std::vector<float>(n_var, 0.0f));
		std::vector<float> XtY(n_var, 0.0f);

		for (size_t k = 0; k < idxs.size(); ++k) {
			const int id = idxs[k];
			std::vector<float> row; row.reserve(n_var);
			row.push_back(1.0f);
			if (use_zwt)   row.push_back((zwt[id] - zwt_min) / std::max(1e-12f, (zwt_max - zwt_min)));
			if (use_dcurv) row.push_back((dcurv[id] - dcurv_min) / std::max(1e-12f, (dcurv_max - dcurv_min)));

			const float y = eq_radius_values[id];
			const float w = out_weights[k];

			for (int j = 0; j < n_var; ++j) {
				XtY[j] += w * row[j] * y;
				for (int l = 0; l < n_var; ++l) XtX[j][l] += w * row[j] * row[l];
			}
		}

		std::vector<std::vector<float>> inv_XtX;
		if (!invert_matrix(XtX, inv_XtX)) return false;

		out_beta.assign(n_var, 0.0f);
		for (int j = 0; j < n_var; ++j)
			for (int l = 0; l < n_var; ++l)
				out_beta[j] += inv_XtX[j][l] * XtY[l];

		return true;
	};

	// ---- 2) Initial fit --------------------------------------------------------
	float zwt_min = 0.0f, zwt_max = 1.0f, dcurv_min = 0.0f, dcurv_max = 1.0f;
	std::vector<float> beta, weights_obs;
	if (!fit_on_subset(valid_indices, use_drift_zwt, use_drift_curv,
		beta, weights_obs, zwt_min, zwt_max, dcurv_min, dcurv_max)) {
		weights_out.assign(N, 0.0f);
		return drift;
	}

	// ---- 3) Geological sign checks --------------------------------------------
	// Expectation: zwt -> beta < 0 ; dcurv -> beta > 0
	bool drift_valid_zwt = use_drift_zwt;
	bool drift_valid_dcurv = use_drift_curv;
	if (use_drift_zwt) {
		const float b_zwt = beta[1];
		if (!(b_zwt < 0.0f)) drift_valid_zwt = false;
	}
	if (use_drift_curv) {
		const int idx = use_drift_zwt ? 2 : 1;
		const float b_dc = beta[idx];
		if (!(b_dc > 0.0f)) drift_valid_dcurv = false;
	}

	auto has_nonzero_slope = [&](const std::vector<float>& b, bool vz, bool vd) -> bool {
		if (vz) { if (std::abs(b[1]) > 1e-6f) return true; }
		if (vd) { const int i2 = vz ? 2 : 1; if (std::abs(b[i2]) > 1e-6f) return true; }
		return false;
	};
	if (!has_nonzero_slope(beta, drift_valid_zwt, drift_valid_dcurv)) {
		weights_out.assign(N, 0.0f);
		return drift;
	}

	// ---- 4) Trimming + optional refit (MAD / Tukey) ---------------------------
	std::vector<float> residuals_cur; residuals_cur.reserve(n_obs);
	for (int id : valid_indices) {
		float yhat = beta[0];
		int bi = 1;
		if (drift_valid_zwt) { const float z01 = (zwt[id] - zwt_min) / std::max(1e-12f, (zwt_max - zwt_min));   yhat += beta[bi++] * z01; }
		if (drift_valid_dcurv) { const float d01 = (dcurv[id] - dcurv_min) / std::max(1e-12f, (dcurv_max - dcurv_min)); yhat += beta[bi++] * d01; }
		residuals_cur.push_back(eq_radius_values[id] - yhat);
	}

	std::vector<float> absr = residuals_cur;
	for (float& v : absr) v = std::abs(v);
	std::nth_element(absr.begin(), absr.begin() + absr.size() / 2, absr.end());
	const float mad = absr[absr.size() / 2];
	const float sigma = std::max(1e-6f, 1.4826f * mad);

	const float C_CUTOFF = 4.685f;
	std::vector<int> survivors; survivors.reserve(n_obs);
	for (int k = 0; k < n_obs; ++k) {
		const float ui = std::abs(residuals_cur[k]) / sigma;
		if (ui <= C_CUTOFF) survivors.push_back(valid_indices[k]);
	}

	const int n_var_cur = 1 + (drift_valid_zwt ? 1 : 0) + (drift_valid_dcurv ? 1 : 0);
	if ((int)survivors.size() >= n_var_cur + 1 && (int)survivors.size() < n_obs) {
		std::vector<float> beta2, weights_obs2;
		float zwt_min2 = zwt_min, zwt_max2 = zwt_max, dcurv_min2 = dcurv_min, dcurv_max2 = dcurv_max;
		if (fit_on_subset(survivors, drift_valid_zwt, drift_valid_dcurv,
			beta2, weights_obs2, zwt_min2, zwt_max2, dcurv_min2, dcurv_max2))
		{
			beta = beta2;
			zwt_min = zwt_min2; zwt_max = zwt_max2;
			dcurv_min = dcurv_min2; dcurv_max = dcurv_max2;
			weights_obs = weights_obs2;
			valid_indices = survivors;
		}
	}

	// ---- 5) Export weights per node (0 for non-observed / trimmed) ------------
	weights_out.assign(N, 0.0f);
	for (size_t k = 0; k < valid_indices.size(); ++k) {
		weights_out[valid_indices[k]] = weights_obs[k];
	}

	// ---- 6) Final drift for all nodes -----------------------------------------
	for (int i = 0; i < N; ++i) {
		float v = beta[0];
		int bi = 1;
		if (drift_valid_zwt) { const float z01 = (zwt[i] - zwt_min) / std::max(1e-12f, (zwt_max - zwt_min));   v += beta[bi++] * z01; }
		if (drift_valid_dcurv) { const float d01 = (dcurv[i] - dcurv_min) / std::max(1e-12f, (dcurv_max - dcurv_min)); v += beta[bi++] * d01; }
		drift[i] = v;
	}

	return drift;
}


void SGS3(
	const KarstNSim::KarsticSkeleton* curve,
	std::vector<float>& simulated_property,
	const std::vector<float>* simulation_distribution,
	const float& global_vario_range,
	const float& global_range_of_neighborhood,
	const float& global_vario_sill,
	const float& global_vario_nugget,
	const std::string& global_vario_model,
	const float& interbranch_vario_range,
	const float& interbranch_range_of_neighborhood,
	const float& interbranch_vario_sill,
	const float& interbranch_vario_nugget,
	const std::string& interbranch_vario_model,
	const float& intrabranch_vario_range,
	const float& intrabranch_range_of_neighborhood,
	const float& intrabranch_vario_sill,
	const float& intrabranch_vario_nugget,
	const std::string& intrabranch_vario_model,
	const int& number_max_of_neighborhood_points,
	const int& nb_points_interbranch,
	const float& proportion_interbranch) {

	// 1) Perform Normal Score Transform of initial distrib AND initial data vector
	std::vector<float> simulated_prop_gauss;
	std::vector<float> simulated_property_copy = simulated_property;
	if (!simulated_property.empty()) {
		simulated_prop_gauss = nst_data_with_nodata(simulated_property, *simulation_distribution, 0., 1.); // gaussianize the data values if any (do not change the no data values though)
	}
	else {
		simulated_prop_gauss.resize(curve->nodes.size());
		std::fill(simulated_prop_gauss.begin(), simulated_prop_gauss.end(), -99999);
	}

	std::vector<float> sim_distrib_gauss(simulation_distribution->size());
	if (!simulation_distribution->empty()) {
		sim_distrib_gauss = nst(*simulation_distribution, 0., 1.);
	}

	// Validate that all conditioning data are inside the support of the initial
	// simulation distribution before launching SGS.
	{
		float sim_min = std::numeric_limits<float>::infinity();
		float sim_max = -std::numeric_limits<float>::infinity();

		for (int i = 0; i < int(simulation_distribution->size()); ++i) {
			const float v = (*simulation_distribution)[i];
			if (std::abs(v - (-99999.0f)) < 1e-12f) {
				continue;
			}
			if (v < sim_min) sim_min = v;
			if (v > sim_max) sim_max = v;
		}

		int n_outside = 0;
		int first_idx = -1;
		float first_val = -99999.0f;

		for (int i = 0; i < int(simulated_property.size()); ++i) {
			const float v = simulated_property[i];

			// Keep only conditioning data
			if (std::abs(v - (-99999.0f)) < 1e-12f) {
				continue;
			}

			if (v < sim_min || v > sim_max) {
				n_outside++;
				if (first_idx < 0) {
					first_idx = i;
					first_val = v;
				}
			}
		}

		if (n_outside > 0) {
			throw std::runtime_error(
				"Sequential Gaussian Simulation aborted: at least one conditioning datum lies outside "
				"the support of the initial simulation distribution. Please change the initial distribution accordingly."
			);
		}
	}

	//save_data(*simulation_distribution, "distrib.txt");
	//save_data(sim_distrib_gauss, "gauss_distrib.txt");

	// 2) Interbranch simulation with interbranch variogram

	// iterate on branches:
	for (int branch_id = 0; branch_id < curve->branch_sizes.size(); branch_id++) {

		// find number of nodes to simulate in given branch
		int nb_nodes_to_simulate = compute_prop_branch(curve->branch_sizes.at(branch_id), nb_points_interbranch, proportion_interbranch);

		std::vector<int> all_branch_nodes;
		// create list of indices of nodes of that branch in nodes, and shuffle that list to get the random order of iteration
		for (int i = 0; i < curve->nodes.size(); i++) {
			if (curve->nodes.at(i).branch_id == branch_id && (simulated_prop_gauss[i] - (-99999)) < 1e-12) { // check that the node is in the same branch and isnt assigned to a value yet (no data value)
				all_branch_nodes.push_back(i);
			}
		}
		std::vector<int> nodes_to_simulate = select_random_elements(all_branch_nodes, nb_nodes_to_simulate);

		for (int i = 0; i < nodes_to_simulate.size(); i++) {

			// determine neighborhood (conditional nodes nearby)

			std::vector<int> neighbors_interbranch = find_neighborhood(nodes_to_simulate.at(i), curve, number_max_of_neighborhood_points, interbranch_range_of_neighborhood, "base", simulated_prop_gauss);

			// if empty neighborhood, simulate directly by sampling from distrib, otherwise, simulate value (kriging + sample from gaussian estimate)
			float val_estimation;
			float var_estimation;

			//kriging_in_point(nodes_to_simulate.at(i), neighbors_interbranch, &sim_distrib_gauss, curve->distance_mat, interbranch_vario_range, interbranch_vario_sill, interbranch_vario_nugget, interbranch_vario_model, simulated_prop_gauss, var_estimation, val_estimation);
			kriging_in_point_on_the_fly(
				nodes_to_simulate.at(i),
				neighbors_interbranch,
				&sim_distrib_gauss,
				curve,
				interbranch_vario_range, interbranch_vario_sill, interbranch_vario_nugget, interbranch_vario_model,
				simulated_prop_gauss,
				var_estimation, val_estimation,
				interbranch_range_of_neighborhood
			);

			if (!neighbors_interbranch.empty()) { // if we had a neighborhood, we simulate value from estimated (kriged) value
				var_estimation = (var_estimation < 0.) ? 0. : var_estimation; // due to numerical uncertainties, var estimation can sometimes be slightly negative.
				float simulated_val = generateNormalRandom(val_estimation, std::sqrt(var_estimation));
				simulated_prop_gauss.at(nodes_to_simulate.at(i)) = simulated_val;
			}
			else { // if the neighborhood was empty, we simply sampled a value from the distribution
				simulated_prop_gauss.at(nodes_to_simulate.at(i)) = val_estimation;
			}
		}
	}

	// 3) Intrabranch simulation with intrabranch variogram

	for (int branch_id = 0; branch_id < curve->branch_sizes.size(); branch_id++) {

		int nb_nodes_to_simulate = compute_prop_branch(curve->branch_sizes.at(branch_id), curve->branch_sizes.at(branch_id), 1.); // compute ALL remaining nodes

		std::vector<int> all_branch_nodes;
		// create list of indices of nodes of that branch in nodes, and shuffle that list to get the random order of iteration
		for (int i = 0; i < curve->nodes.size(); i++) {
			if (curve->nodes.at(i).branch_id == branch_id && (simulated_prop_gauss[i] - (-99999)) < 1e-12) { // check that the node is in the same branch and isnt assigned to a value yet (no data value)
				all_branch_nodes.push_back(i);
			}
		}
		std::vector<int> nodes_to_simulate = select_random_elements(all_branch_nodes, nb_nodes_to_simulate);

		for (int i = 0; i < nodes_to_simulate.size(); i++) {
			std::vector<int> neighbors_intrabranch = find_neighborhood(nodes_to_simulate.at(i), curve, number_max_of_neighborhood_points, intrabranch_range_of_neighborhood, "branch", simulated_prop_gauss);

			float val_estimation;
			float var_estimation;

			//kriging_in_point(nodes_to_simulate.at(i), neighbors_intrabranch, &sim_distrib_gauss, curve->distance_mat, intrabranch_vario_range, intrabranch_vario_sill, intrabranch_vario_nugget, intrabranch_vario_model, simulated_prop_gauss, var_estimation, val_estimation);
			kriging_in_point_on_the_fly(
				nodes_to_simulate.at(i),
				neighbors_intrabranch,
				&sim_distrib_gauss,
				curve,
				intrabranch_vario_range, intrabranch_vario_sill, intrabranch_vario_nugget, intrabranch_vario_model,
				simulated_prop_gauss,
				var_estimation, val_estimation,
				intrabranch_range_of_neighborhood
			);


			if (!neighbors_intrabranch.empty()) {
				var_estimation = (var_estimation < 0.) ? 0. : var_estimation;// due to numerical uncertainties, var estimation can sometimes be slightly negative.
				float simulated_val = generateNormalRandom(val_estimation, std::sqrt(var_estimation));
				simulated_prop_gauss.at(nodes_to_simulate.at(i)) = simulated_val;
			}
			else {
				simulated_prop_gauss.at(nodes_to_simulate.at(i)) = val_estimation;
			}
		}
	}

	// 4) Intersections simulation with global variogram

	std::vector<int> all_intersection_nodes;
	for (int i = 0; i < curve->nodes.size(); i++) {
		if (curve->nodes.at(i).branch_id == -1 && (simulated_prop_gauss[i] - (-99999)) < 1e-12) // if intersection and value not already assigned
			all_intersection_nodes.push_back(i);
	}
	std::vector<int> nodes_to_simulate_global = select_random_elements(all_intersection_nodes, int(all_intersection_nodes.size())); // we compute ALL intersections, no exceptions

	for (int i = 0; i < nodes_to_simulate_global.size(); i++) {

		std::vector<int> neighbors_global = find_neighborhood(nodes_to_simulate_global.at(i), curve, number_max_of_neighborhood_points, global_range_of_neighborhood, "base", simulated_prop_gauss);

		float val_estimation_global;
		float var_estimation_global;

		//kriging_in_point(nodes_to_simulate_global.at(i), neighbors_global, &sim_distrib_gauss, curve->distance_mat, global_vario_range, global_vario_sill, global_vario_nugget, global_vario_model, simulated_prop_gauss, var_estimation_global, val_estimation_global);
		kriging_in_point_on_the_fly(
			nodes_to_simulate_global.at(i),
			neighbors_global,
			&sim_distrib_gauss,
			curve,
			global_vario_range, global_vario_sill, global_vario_nugget, global_vario_model,
			simulated_prop_gauss,
			var_estimation_global, val_estimation_global,
			global_range_of_neighborhood
		);

		if (!neighbors_global.empty()) {
			var_estimation_global = (var_estimation_global < 0.) ? 0. : var_estimation_global;// due to numerical uncertainties, var estimation can sometimes be slightly negative.
			float simulated_val_global = generateNormalRandom(val_estimation_global, std::sqrt(var_estimation_global));
			simulated_prop_gauss.at(nodes_to_simulate_global.at(i)) = simulated_val_global;
		}
		else {
			simulated_prop_gauss.at(nodes_to_simulate_global.at(i)) = val_estimation_global;
		}
	}

	// 5) Back-transform of generated distribution.

	// if we had data at the beginning, use it as basis for back transformation
	//if (!simulated_property.empty()) {
	//	simulated_property = back_transform(simulated_prop_gauss, simulated_property_copy);
	//}
	//else { // else, use the initial distribution instead
	simulated_property = back_transform(simulated_prop_gauss, *simulation_distribution);
	//}
}

int compute_prop_branch(const int& total_nb_nodes_branch, const int& nb_points_interbranch, const float& proportion_interbranch) {

	float true_proportion_interbranch;
	int needed_value_points_for_this_branch;

	if (proportion_interbranch > 1.) {
		true_proportion_interbranch = 1.;
	}
	else {
		true_proportion_interbranch = proportion_interbranch;
	}

	int proportion_value = round(total_nb_nodes_branch * true_proportion_interbranch);

	if (proportion_value <= nb_points_interbranch) {
		needed_value_points_for_this_branch = proportion_value;
	}
	else {
		needed_value_points_for_this_branch = nb_points_interbranch;
	}

	return needed_value_points_for_this_branch;
}

// --- SGS3_with_external_drift ---
// Sequential Gaussian Simulation with an external drift field.
// Workflow:
//   1. Compute external drift field m(x) using weighted OLS regression.
//   2. If drift is degenerate (~0 everywhere), fallback to plain SGS3.
//   3. Compute residuals on observed nodes: Re(x) - m(x).
//   4. Simulate residuals with SGS3, using a zero-mean distribution
//      to avoid adding a constant bias at recomposition.
//   5. Recompose the final property: simulated_property = residuals + drift.
void SGS3_with_external_drift(
	const KarstNSim::KarsticSkeleton* curve,
	const std::vector<Vector3>& springs_xyz,
	std::vector<float>& simulated_property,
	const std::vector<float>* simulation_distribution,
	const float& global_vario_range,
	const float& global_range_of_neighborhood,
	const float& global_vario_sill,
	const float& global_vario_nugget,
	const std::string& global_vario_model,
	const float& interbranch_vario_range,
	const float& interbranch_range_of_neighborhood,
	const float& interbranch_vario_sill,
	const float& interbranch_vario_nugget,
	const std::string& interbranch_vario_model,
	const float& intrabranch_vario_range,
	const float& intrabranch_range_of_neighborhood,
	const float& intrabranch_vario_sill,
	const float& intrabranch_vario_nugget,
	const std::string& intrabranch_vario_model,
	const int& number_max_of_neighborhood_points,
	const int& nb_points_interbranch,
	const float& proportion_interbranch,
	const float& z_phreatic,
	const bool& use_drift_zwt,
	const bool& use_drift_curv,
	std::vector<float>& drift_output,
	std::vector<float>& weights_output)
{
	// === 1) Compute external drift field m(x) ===
	std::vector<float> drift = compute_external_drift(
		curve,
		springs_xyz,
		z_phreatic,
		use_drift_zwt,
		use_drift_curv,
		simulated_property,
		weights_output
	);
	drift_output = drift;
	// === 2) If drift is degenerate (all values ~0), fallback to plain SGS3 ===
	bool drift_disabled = std::all_of(drift.begin(), drift.end(),
		[](float v) { return std::abs(v) < 1e-8f; });
	if (drift_disabled) {
			SGS3(
				curve,
				simulated_property,
				simulation_distribution,
				global_vario_range,
				global_range_of_neighborhood,
				global_vario_sill,
				global_vario_nugget,
				global_vario_model,
				interbranch_vario_range,
				interbranch_range_of_neighborhood,
				interbranch_vario_sill,
				interbranch_vario_nugget,
				interbranch_vario_model,
				intrabranch_vario_range,
				intrabranch_range_of_neighborhood,
				intrabranch_vario_sill,
				intrabranch_vario_nugget,
				intrabranch_vario_model,
				number_max_of_neighborhood_points,
				nb_points_interbranch,
				proportion_interbranch
			);
		return; // no recomposition residual + drift
	}
	// === 3) Compute residuals on observed nodes: Re(x) - m(x) ===
	std::vector<float> residuals(simulated_property.size(), -99999.0f);
	for (size_t i = 0; i < simulated_property.size(); ++i) {
		if (std::abs(simulated_property[i] - (-99999.0f)) > 1e-12f) {
			residuals[i] = simulated_property[i] - drift[i];
		}
	}
	// === 4) Simulate residuals with a zero-centered distribution ============
	// Rationale:
	//  - If the provided simulation distribution has a non-zero mean,
	//    the simulated residuals would inherit that bias and add it to m(x),
	//    creating an artificial offset in the final property.
	//  - To avoid this, we re-center the distribution to zero mean before
	//    passing it to SGS3. This ensures residuals have mean ~0.
	std::vector<float> residual_sim_distribution;
	if (simulation_distribution && !simulation_distribution->empty()) {
		residual_sim_distribution = *simulation_distribution;
		double mu = std::accumulate(residual_sim_distribution.begin(),
			residual_sim_distribution.end(), 0.0)
			/ static_cast<double>(residual_sim_distribution.size());
		for (float& v : residual_sim_distribution) v -= static_cast<float>(mu);

		// Safety: if everything cancels to ~0, keep a tiny symmetric spread
		bool degenerate = true;
		for (float v : residual_sim_distribution) {
			if (std::abs(v) > 1e-12f) { degenerate = false; break; }
		}
		if (degenerate) residual_sim_distribution = { -1.f, 0.f, 1.f };
	}
	else {
		// No distribution provided -> create a small zero-mean placeholder
		residual_sim_distribution = { -1.f, 0.f, 1.f };
	}

	SGS3(
		curve,
		residuals,
		&residual_sim_distribution, // zero-mean residual distribution
		global_vario_range,
		global_range_of_neighborhood,
		global_vario_sill,
		global_vario_nugget,
		global_vario_model,
		interbranch_vario_range,
		interbranch_range_of_neighborhood,
		interbranch_vario_sill,
		interbranch_vario_nugget,
		interbranch_vario_model,
		intrabranch_vario_range,
		intrabranch_range_of_neighborhood,
		intrabranch_vario_sill,
		intrabranch_vario_nugget,
		intrabranch_vario_model,
		number_max_of_neighborhood_points,
		nb_points_interbranch,
		proportion_interbranch
	);

	// === 5) Recompose final property: simulated_property = residuals + drift ===
	simulated_property.resize(residuals.size());
	for (size_t i = 0; i < residuals.size(); ++i) {
		if (std::abs(residuals[i] - (-99999.0f)) > 1e-12) {
			simulated_property[i] = residuals[i] + drift[i];
		}
		else {
			simulated_property[i] = -99999.0f; // untouched NDV
		}
	}
}


