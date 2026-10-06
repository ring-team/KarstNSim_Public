/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr for new methods + modifications to original methods
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

The base SGS3 algorithm is rewritten (except the variogram_value fonction, written by G. Rongier in 2015) from the pseudo-code algorithm (Algorithm 2) in the paper of Frantz et al., 2021 "Analysis and stochastic simulation of geometrical properties of conduits in karstic networks".
This work was performed in the frame of the RING project at Université de Lorraine.

***************************************************************/

#include "KarstNSim/geostats.h"
#include <array>
#include <cstdint>
#include <unordered_map>


namespace {
	// Variogram-parameter convention used by SGS:
	// false -> sill and nugget are supplied in the simulated-property space and
	//          are converted internally to the same normal-score space used here.
	// true  -> sill and nugget are already expressed in normal-score space and
	//          are used directly
	constexpr bool K_VARIOGRAM_PARAMETERS_ARE_ALREADY_GAUSSIAN = true;

	// Enables the equivalent-radius upper threshold during external-drift regression:
	// true  -> observations with an equivalent radius greater than
	//          K_RADIUS_MAX_FOR_REGRESSION are excluded from the regression.
	// false -> all valid equivalent-radius observations are used, regardless of radius.
	constexpr bool K_EXT_DRIFT_ENABLE_RADIUS_CAP = false;

	// Enables robust outlier rejection during external-drift regression:
	// true  -> observations with excessive MAD-standardized residuals may be rejected,
	//          followed by a regression refit using the retained observations.
	// false -> no MAD-based observation trimming or subsequent refit is performed.
	constexpr bool K_EXT_DRIFT_ENABLE_MAD_TRIMMING = true;

	// Enables compact diagnostic logging for the external-drift SGS workflow:
	// true  -> reports regression filtering, fitted coefficients, robust trimming,
	//          predictor ranges, final drift range, and residual-distribution scaling.
	// false -> no additional SGS diagnostic output is produced.
	constexpr bool K_EXT_DRIFT_ENABLE_DIAGNOSTIC_LOGS = false;

	// Maximum equivalent radius accepted in the external-drift regression when
	// K_EXT_DRIFT_ENABLE_RADIUS_CAP is true. The value uses the same length unit
	// as the equivalent-radius conditioning data.
	constexpr float K_RADIUS_MAX_FOR_REGRESSION = 2.1f;

	// Parameters used for variogram conversion to normal-score space, if
	// K_VARIOGRAM_PARAMETERS_ARE_ALREADY_GAUSSIAN is false.
	constexpr int K_TRANS_GAUSSIAN_HERMITE_TERMS = 24;
	constexpr int K_TRANS_GAUSSIAN_LOOKUP_SIZE = 4097;

	/**
	 * @brief Approximates the inverse standard-normal cumulative distribution.
	 *
	 * The rational approximation is deterministic and sufficiently accurate for
	 * constructing the empirical anamorphosis at probability midpoints.
	 */
	double inverse_standard_normal_cdf(const double probability)
	{
		if (!(probability > 0.0 && probability < 1.0)) {
			throw std::invalid_argument(
				"Normal-score probabilities must lie strictly between zero and one.");
		}

		const double a1 = -3.969683028665376e+01;
		const double a2 = 2.209460984245205e+02;
		const double a3 = -2.759285104469687e+02;
		const double a4 = 1.383577518672690e+02;
		const double a5 = -3.066479806614716e+01;
		const double a6 = 2.506628277459239e+00;
		const double b1 = -5.447609879822406e+01;
		const double b2 = 1.615858368580409e+02;
		const double b3 = -1.556989798598866e+02;
		const double b4 = 6.680131188771972e+01;
		const double b5 = -1.328068155288572e+01;
		const double c1 = -7.784894002430293e-03;
		const double c2 = -3.223964580411365e-01;
		const double c3 = -2.400758277161838e+00;
		const double c4 = -2.549732539343734e+00;
		const double c5 = 4.374664141464968e+00;
		const double c6 = 2.938163982698783e+00;
		const double d1 = 7.784695709041462e-03;
		const double d2 = 3.224671290700398e-01;
		const double d3 = 2.445134137142996e+00;
		const double d4 = 3.754408661907416e+00;
		const double lower_tail = 0.02425;
		const double upper_tail = 1.0 - lower_tail;

		if (probability < lower_tail) {
			const double q = std::sqrt(-2.0 * std::log(probability));
			return (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
				((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0);
		}
		if (probability > upper_tail) {
			const double q = std::sqrt(-2.0 * std::log(1.0 - probability));
			return -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
				((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0);
		}

		const double q = probability - 0.5;
		const double r = q * q;
		return (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q /
			(((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1.0);
	}

	/**
	 * @brief Converts property-space semivariances to normal-score semivariances.
	 *
	 * The empirical anamorphosis is expanded on normalized Hermite polynomials.
	 * For a Gaussian pair with correlation rho, the expansion gives the
	 * property-space covariance as a positive power series in rho. Inverting
	 * that monotone relation yields the normal-score semivariance 1-rho.
	 */
	class NormalScoreVariogramConverter {
	public:
		explicit NormalScoreVariogramConverter(const std::vector<float>& distribution)
		{
			std::vector<double> values;
			values.reserve(distribution.size());
			for (const float value : distribution) {
				if (!std::isfinite(value)) {
					throw std::invalid_argument(
						"The simulation distribution contains a non-finite value.");
				}
				values.push_back(static_cast<double>(value));
			}
			if (values.size() < 2) {
				throw std::invalid_argument(
					"At least two finite values are required to convert variograms to normal-score space.");
			}

			std::sort(values.begin(), values.end());
			const double count = static_cast<double>(values.size());
			const double mean = std::accumulate(values.begin(), values.end(), 0.0) / count;
			property_variance_ = 0.0;
			for (const double value : values) {
				const double centered = value - mean;
				property_variance_ += centered * centered;
			}
			property_variance_ /= count;
			if (!std::isfinite(property_variance_) ||
				property_variance_ <= std::numeric_limits<double>::epsilon()) {
				throw std::invalid_argument(
					"The simulation distribution has zero or non-finite variance.");
			}

			const int term_count = std::min(
				K_TRANS_GAUSSIAN_HERMITE_TERMS,
				std::max(1, static_cast<int>(values.size()) - 1));
			std::vector<double> coefficients(term_count + 1, 0.0);

			for (size_t sample = 0; sample < values.size(); ++sample) {
				const double probability =
					(static_cast<double>(sample) + 0.5) / count;
				const double gaussian_score = inverse_standard_normal_cdf(probability);
				const double centered_value = values[sample] - mean;

				double psi_previous = 1.0;
				double psi_current = gaussian_score;
				coefficients[1] += centered_value * psi_current / count;

				for (int order = 2; order <= term_count; ++order) {
					const double psi_next =
						(gaussian_score * psi_current -
							std::sqrt(static_cast<double>(order - 1)) * psi_previous) /
						std::sqrt(static_cast<double>(order));
					coefficients[order] += centered_value * psi_next / count;
					psi_previous = psi_current;
					psi_current = psi_next;
				}
			}

			covariance_terms_.assign(term_count + 1, 0.0);
			double represented_variance = 0.0;
			for (int order = 1; order <= term_count; ++order) {
				covariance_terms_[order] = coefficients[order] * coefficients[order];
				represented_variance += covariance_terms_[order];
			}
			if (represented_variance <= std::numeric_limits<double>::epsilon()) {
				covariance_terms_.assign(2, 0.0);
				covariance_terms_[1] = property_variance_;
			}
			else {
				const double variance_scale = property_variance_ / represented_variance;
				for (size_t order = 1; order < covariance_terms_.size(); ++order) {
					covariance_terms_[order] *= variance_scale;
				}
			}

			build_lookup_table();
		}

		double gaussian_semivariance(const double property_semivariance) const
		{
			if (property_semivariance <= 0.0) return 0.0;
			if (property_semivariance >= property_variance_) return 1.0;

			const double position = property_semivariance / property_variance_ *
				static_cast<double>(gaussian_semivariance_lookup_.size() - 1);
			const size_t lower = static_cast<size_t>(position);
			const size_t upper = std::min(
				lower + 1, gaussian_semivariance_lookup_.size() - 1);
			const double fraction = position - static_cast<double>(lower);
			return gaussian_semivariance_lookup_[lower] * (1.0 - fraction) +
				gaussian_semivariance_lookup_[upper] * fraction;
		}

	private:
		double covariance_for_correlation(const double correlation) const
		{
			double covariance = 0.0;
			for (size_t order = covariance_terms_.size() - 1; order > 0; --order) {
				covariance = (covariance + covariance_terms_[order]) * correlation;
			}
			return covariance;
		}

		void build_lookup_table()
		{
			gaussian_semivariance_lookup_.assign(
				K_TRANS_GAUSSIAN_LOOKUP_SIZE, 0.0);
			gaussian_semivariance_lookup_.back() = 1.0;

			for (int index = 1; index < K_TRANS_GAUSSIAN_LOOKUP_SIZE - 1; ++index) {
				const double target_semivariance = property_variance_ *
					static_cast<double>(index) /
					static_cast<double>(K_TRANS_GAUSSIAN_LOOKUP_SIZE - 1);
				double lower_correlation = 0.0;
				double upper_correlation = 1.0;

				for (int iteration = 0; iteration < 56; ++iteration) {
					const double correlation =
						0.5 * (lower_correlation + upper_correlation);
					const double candidate_semivariance =
						property_variance_ - covariance_for_correlation(correlation);
					if (candidate_semivariance > target_semivariance) {
						lower_correlation = correlation;
					}
					else {
						upper_correlation = correlation;
					}
				}

				gaussian_semivariance_lookup_[index] =
					1.0 - 0.5 * (lower_correlation + upper_correlation);
			}
		}

		double property_variance_ = 0.0;
		std::vector<double> covariance_terms_;
		std::vector<double> gaussian_semivariance_lookup_;
	};

	thread_local const NormalScoreVariogramConverter* active_variogram_converter = nullptr;

	class ScopedVariogramConverter {
	public:
		explicit ScopedVariogramConverter(const NormalScoreVariogramConverter* converter)
			: previous_(active_variogram_converter)
		{
			active_variogram_converter = converter;
		}

		~ScopedVariogramConverter()
		{
			active_variogram_converter = previous_;
		}

		ScopedVariogramConverter(const ScopedVariogramConverter&) = delete;
		ScopedVariogramConverter& operator=(const ScopedVariogramConverter&) = delete;

	private:
		const NormalScoreVariogramConverter* previous_;
	};

	double active_zero_lag_variance(const double supplied_sill)
	{
		return active_variogram_converter == nullptr ? supplied_sill : 1.0;
	}

	/**
	 * @brief Reusable scratch storage for truncated Dijkstra searches.
	 *
	 * Generation counters avoid clearing and reallocating arrays with one entry per
	 * skeleton node for every local geostatistical search. Only entries touched by
	 * the current search are logically initialized. This is strictly equivalent to
	 * rebuilding distance, settled, and target arrays at each call.
	 */
	struct GeostatsDijkstraWorkspace {
		std::vector<float> distance;
		std::vector<std::uint32_t> distance_generation;
		std::vector<std::uint32_t> settled_generation;
		std::vector<std::uint32_t> target_generation;
		std::uint32_t generation = 0;

		/**
		 * @brief Starts a new logical search and resizes storage when required.
		 * @param node_count Number of skeleton nodes.
		 */
		void begin(const std::size_t node_count)
		{
			if (distance.size() != node_count) {
				distance.assign(node_count, 0.0f);
				distance_generation.assign(node_count, 0u);
				settled_generation.assign(node_count, 0u);
				target_generation.assign(node_count, 0u);
				generation = 1u;
				return;
			}

			++generation;
			if (generation == 0u) {
				std::fill(distance_generation.begin(), distance_generation.end(), 0u);
				std::fill(settled_generation.begin(), settled_generation.end(), 0u);
				std::fill(target_generation.begin(), target_generation.end(), 0u);
				generation = 1u;
			}
		}

		/**
		 * @brief Returns the current-search distance or infinity when untouched.
		 * @param node Skeleton node index.
		 * @return Current tentative distance.
		 */
		float get_distance(const int node) const
		{
			return distance_generation[static_cast<std::size_t>(node)] == generation
				? distance[static_cast<std::size_t>(node)]
				: std::numeric_limits<float>::infinity();
		}

		/**
		 * @brief Stores a tentative distance for the current search.
		 * @param node Skeleton node index.
		 * @param value Tentative shortest-path distance.
		 */
		void set_distance(const int node, const float value)
		{
			const std::size_t index = static_cast<std::size_t>(node);
			distance[index] = value;
			distance_generation[index] = generation;
		}

		/**
		 * @brief Tests whether a node has already been settled in the current search.
		 * @param node Skeleton node index.
		 * @return True when the node is settled.
		 */
		bool is_settled(const int node) const
		{
			return settled_generation[static_cast<std::size_t>(node)] == generation;
		}

		/**
		 * @brief Marks a node as settled in the current search.
		 * @param node Skeleton node index.
		 */
		void mark_settled(const int node)
		{
			settled_generation[static_cast<std::size_t>(node)] = generation;
		}

		/**
		 * @brief Marks a node as a target in the current search.
		 * @param node Skeleton node index.
		 */
		void mark_target(const int node)
		{
			target_generation[static_cast<std::size_t>(node)] = generation;
		}

		/**
		 * @brief Tests whether a node is a target in the current search.
		 * @param node Skeleton node index.
		 * @return True when the node is a target.
		 */
		bool is_target(const int node) const
		{
			return target_generation[static_cast<std::size_t>(node)] == generation;
		}
	};

	thread_local GeostatsDijkstraWorkspace GEOSTATS_DIJKSTRA_WORKSPACE;

}

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

	auto& workspace = GEOSTATS_DIJKSTRA_WORKSPACE;
	workspace.begin(static_cast<std::size_t>(N));

	auto edge_length = [&](int u, int v) -> float {
		// Geometric edge weight; same logic as before (Euclidean length).
		const Vector3& pu = curve->nodes[u].p;
		const Vector3& pv = curve->nodes[v].p;
		return KarstNSim::magnitude(pu - pv);
	};

	std::priority_queue<NodeKey, std::vector<NodeKey>, std::greater<NodeKey>> pq;
	workspace.set_distance(current_node_index, 0.f);
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
		if (workspace.is_settled(u)) continue;
		workspace.mark_settled(u);

		// Range truncation: in Dijkstra, if the minimum key exceeds R, remaining keys are >= du > R.
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
			if (workspace.is_settled(v)) continue;

			const float w = edge_length(u, v);
			if (!std::isfinite(w) || w < 0.f) continue;

			const float alt = du + w;
			if (!std::isfinite(alt) || alt > range_of_neighborhood) continue;

			if (alt < workspace.get_distance(v)) {
				workspace.set_distance(v, alt);
				pq.push({ alt, v });
			}
		}
	}

	// Sort by metric distance and truncate to K nearest within range.
	// The original full sort is deliberately retained because its tie behavior is
	// part of the existing simulation path and therefore must remain unchanged.
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

	if (vario_value >= 0.0f && active_variogram_converter != nullptr) {
		vario_value = static_cast<float>(
			active_variogram_converter->gaussian_semivariance(vario_value));
	}
	return vario_value;
}

// Function to calculate the inverse of a square matrix using Gaussian elimination
bool invert_matrix(const std::vector<std::vector<float>>& input, std::vector<std::vector<float>>& output) {
	if (input.empty()) {
		std::cerr << "Cannot invert an empty matrix." << std::endl;
		return false;
	}

	// Check if the input matrix is square
	int size = int(input.size());
	if (size != int(input[0].size())) {
		std::cerr << "Input matrix is not square." << std::endl;
		return false;
	}
	for (const auto& row : input) {
		if (int(row.size()) != size) {
			std::cerr << "Input matrix is not square." << std::endl;
			return false;
		}
	}

	float matrix_scale = 0.0f;
	for (const auto& row : input) {
		for (float value : row) {
			if (!std::isfinite(value)) {
				std::cerr << "Input matrix contains a non-finite value." << std::endl;
				return false;
			}
			matrix_scale = std::max(matrix_scale, std::abs(value));
		}
	}
	if (matrix_scale == 0.0f) {
		std::cerr << "Matrix is singular: all coefficients are zero." << std::endl;
		return false;
	}
	const float pivot_tolerance =
		std::numeric_limits<float>::epsilon() *
		static_cast<float>(std::max(1, size)) * matrix_scale;

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

		// Detect numerical singularity using a scale-aware threshold rather than
		// testing only for an exactly zero pivot.
		if (!std::isfinite(A[i][i]) || std::abs(A[i][i]) <= pivot_tolerance) {
			std::cerr << "Matrix is numerically singular at pivot " << i
				<< " (|pivot|=" << std::abs(A[i][i])
				<< ", tolerance=" << pivot_tolerance << ")." << std::endl;
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

namespace {
	/**
	 * @brief Solves a symmetric positive-definite linear system by Cholesky decomposition.
	 *
	 * The matrix is factorized as A = L L^T, followed by forward and backward
	 * substitutions. All linear-algebra operations are performed in double
	 * precision. A scale-aware pivot tolerance is used to detect matrices that
	 * are singular or numerically non-positive-definite.
	 *
	 * @param matrix Symmetric positive-definite coefficient matrix.
	 * @param rhs Right-hand-side vector.
	 * @param solution Solution vector, overwritten on success.
	 * @param failure_reason Concise diagnostic filled on failure.
	 * @return true when the system was solved successfully; false otherwise.
	 */
	bool solve_spd_cholesky(
		const std::vector<std::vector<double>>& matrix,
		const std::vector<double>& rhs,
		std::vector<double>& solution,
		std::string& failure_reason)
	{
		const int size = static_cast<int>(matrix.size());
		if (size == 0 || static_cast<int>(rhs.size()) != size) {
			failure_reason = "empty system or incompatible right-hand-side size";
			return false;
		}

		double matrix_scale = 0.0;
		for (const auto& row : matrix) {
			if (static_cast<int>(row.size()) != size) {
				failure_reason = "coefficient matrix is not square";
				return false;
			}
			for (double value : row) {
				if (!std::isfinite(value)) {
					failure_reason = "coefficient matrix contains a non-finite value";
					return false;
				}
				matrix_scale = std::max(matrix_scale, std::abs(value));
			}
		}
		for (double value : rhs) {
			if (!std::isfinite(value)) {
				failure_reason = "right-hand side contains a non-finite value";
				return false;
			}
		}
		if (matrix_scale == 0.0) {
			failure_reason = "covariance matrix is identically zero";
			return false;
		}

		const double pivot_tolerance =
			std::numeric_limits<double>::epsilon() *
			static_cast<double>(std::max(1, size)) * matrix_scale;
		std::vector<std::vector<double>> lower(
			size, std::vector<double>(size, 0.0));

		for (int i = 0; i < size; ++i) {
			for (int j = 0; j <= i; ++j) {
				double value = matrix[i][j];
				for (int k = 0; k < j; ++k) {
					value -= lower[i][k] * lower[j][k];
				}

				if (i == j) {
					if (!std::isfinite(value) || value <= pivot_tolerance) {
						std::ostringstream diagnostic;
						diagnostic << "Cholesky factorization failed at pivot " << i
							<< " (value=" << value
							<< ", tolerance=" << pivot_tolerance << ")";
						failure_reason = diagnostic.str();
						return false;
					}
					lower[i][j] = std::sqrt(value);
				}
				else {
					lower[i][j] = value / lower[j][j];
				}
			}
		}

		std::vector<double> intermediate(size, 0.0);
		for (int i = 0; i < size; ++i) {
			double value = rhs[i];
			for (int j = 0; j < i; ++j) {
				value -= lower[i][j] * intermediate[j];
			}
			intermediate[i] = value / lower[i][i];
		}

		solution.assign(size, 0.0);
		for (int i = size - 1; i >= 0; --i) {
			double value = intermediate[i];
			for (int j = i + 1; j < size; ++j) {
				value -= lower[j][i] * solution[j];
			}
			solution[i] = value / lower[i][i];
			if (!std::isfinite(solution[i])) {
				failure_reason = "linear solve produced a non-finite kriging weight";
				return false;
			}
		}

		return true;
	}
}


namespace {
	/**
	 * @brief Accelerates exact leave-one-out weighted regressions used by drift trimming.
	 *
	 * The redundancy kernel depends only on a small local neighborhood: at most
	 * 17 samples in the two-predictor case and at most 2*K+1 samples in the
	 * one-predictor case, with K <= 8. Removing one observation therefore changes
	 * only the redundancy weights of observations whose local neighborhood contains
	 * the removed sample, provided that predictor normalization bounds are unchanged.
	 *
	 * This cache stores those local neighborhoods once. For each eligible leave-one-
	 * out fit it recomputes exactly the affected raw weights, then performs the same
	 * min/max weight rescaling, class balancing, weighted normal-equation assembly,
	 * and matrix inversion in the same observation order as the original code.
	 * Exclusions that change a predictor minimum/maximum, or a neighborhood-size
	 * threshold, are intentionally rejected by this fast path and must use the
	 * original full-fit fallback.
	 */
	class FastLeaveOneOutRegressionCache {
	public:
		/**
		 * @brief Builds the reusable local-neighborhood representation.
		 * @param indices Regression observation node indices in fitting order.
		 * @param zwt Predictor values for vertical distance above phreatic level.
		 * @param dcurv Predictor values for upstream curvilinear length.
		 * @param response Equivalent-radius response values indexed by skeleton node.
		 * @param is_spring_node Precomputed spring-class flag indexed by skeleton node.
		 * @param use_zwt Whether zwt is active in the current regression.
		 * @param use_dcurv Whether dcurv is active in the current regression.
		 * @param zwt_min Normalization minimum used by the full current fit.
		 * @param zwt_max Normalization maximum used by the full current fit.
		 * @param dcurv_min Normalization minimum used by the full current fit.
		 * @param dcurv_max Normalization maximum used by the full current fit.
		 */
		FastLeaveOneOutRegressionCache(
			const std::vector<int>& indices,
			const std::vector<float>& zwt,
			const std::vector<float>& dcurv,
			const std::vector<float>& response,
			const std::vector<std::uint8_t>& is_spring_node,
			const bool use_zwt,
			const bool use_dcurv,
			const float zwt_min,
			const float zwt_max,
			const float dcurv_min,
			const float dcurv_max)
			: indices_(indices),
			zwt_(zwt),
			dcurv_(dcurv),
			response_(response),
			is_spring_node_(is_spring_node),
			use_zwt_(use_zwt),
			use_dcurv_(use_dcurv),
			zwt_min_(zwt_min),
			zwt_max_(zwt_max),
			dcurv_min_(dcurv_min),
			dcurv_max_(dcurv_max)
		{
			initialize();
		}

		/**
		 * @brief Tests whether one exclusion can use the cached exact fast path.
		 * @param excluded_position Position of the removed observation in `indices`.
		 * @return True when normalization and neighborhood-size rules are unchanged.
		 */
		bool can_fit_excluding(const int excluded_position) const
		{
			if (!enabled_ || excluded_position < 0 || excluded_position >= observation_count_) {
				return false;
			}

			const int node_id = indices_[static_cast<std::size_t>(excluded_position)];
			if (use_zwt_ && exclusion_changes_range(
				zwt_[static_cast<std::size_t>(node_id)],
				zwt_actual_min_, zwt_actual_max_, zwt_min_count_, zwt_max_count_)) {
				return false;
			}
			if (use_dcurv_ && exclusion_changes_range(
				dcurv_[static_cast<std::size_t>(node_id)],
				dcurv_actual_min_, dcurv_actual_max_, dcurv_min_count_, dcurv_max_count_)) {
				return false;
			}
			return true;
		}

		/**
		 * @brief Fits the weighted regression after removing one observation.
		 * @param excluded_position Position of the removed observation in `indices`.
		 * @param beta Regression coefficients overwritten on success.
		 * @return True when the cached fit is valid and the normal matrix is invertible.
		 */
		bool fit_excluding(const int excluded_position, std::vector<float>& beta)
		{
			if (!can_fit_excluding(excluded_position)) return false;
			begin_override_generation();

			const std::vector<int>& affected =
				affected_by_exclusion_[static_cast<std::size_t>(excluded_position)];
			for (int position : affected) {
				if (position == excluded_position) continue;
				override_raw_[static_cast<std::size_t>(position)] =
					recompute_raw_weight(position, excluded_position);
				override_generation_[static_cast<std::size_t>(position)] = override_generation_id_;
			}

			float wmin = std::numeric_limits<float>::infinity();
			float wmax = -std::numeric_limits<float>::infinity();

			for (const auto& item : raw_weight_order_) {
				const int position = item.second;
				if (position == excluded_position || has_override(position)) continue;
				wmin = item.first;
				break;
			}
			for (auto it = raw_weight_order_.rbegin(); it != raw_weight_order_.rend(); ++it) {
				const int position = it->second;
				if (position == excluded_position || has_override(position)) continue;
				wmax = it->first;
				break;
			}
			for (int position : affected) {
				if (position == excluded_position || !has_override(position)) continue;
				const float value = override_raw_[static_cast<std::size_t>(position)];
				wmin = std::min(wmin, value);
				wmax = std::max(wmax, value);
			}

			if (!std::isfinite(wmin) || !std::isfinite(wmax)) return false;

			float Wspr = 0.0f;
			float Woth = 0.0f;
			for (int position = 0; position < observation_count_; ++position) {
				if (position == excluded_position) continue;
				const float raw = raw_weight_at(position);
				float scaled = 1.0f;
				if (wmax > wmin) {
					scaled = 0.2f + 0.8f * (raw - wmin) / (wmax - wmin);
				}
				scaled_weight_workspace_[static_cast<std::size_t>(position)] = scaled;

				const int node_id = indices_[static_cast<std::size_t>(position)];
				if (is_spring_node_[static_cast<std::size_t>(node_id)] != 0u) Wspr += scaled;
				else Woth += scaled;
			}

			const float Wtot = Wspr + Woth;
			if (Wspr > 0.0f && Woth > 0.0f) {
				const float target = 0.5f * Wtot;
				const float fs = target / Wspr;
				const float fo = target / Woth;
				for (int position = 0; position < observation_count_; ++position) {
					if (position == excluded_position) continue;
					const int node_id = indices_[static_cast<std::size_t>(position)];
					scaled_weight_workspace_[static_cast<std::size_t>(position)] *=
						(is_spring_node_[static_cast<std::size_t>(node_id)] != 0u ? fs : fo);
				}
			}

			const int n_var = 1 + (use_zwt_ ? 1 : 0) + (use_dcurv_ ? 1 : 0);
			std::vector<std::vector<float>> XtX(
				static_cast<std::size_t>(n_var),
				std::vector<float>(static_cast<std::size_t>(n_var), 0.0f));
			std::vector<float> XtY(static_cast<std::size_t>(n_var), 0.0f);

			for (int position = 0; position < observation_count_; ++position) {
				if (position == excluded_position) continue;
				const int node_id = indices_[static_cast<std::size_t>(position)];

				std::vector<float> row;
				row.reserve(static_cast<std::size_t>(n_var));
				row.push_back(1.0f);
				if (use_zwt_) {
					row.push_back(
						(zwt_[static_cast<std::size_t>(node_id)] - zwt_min_) /
						std::max(1e-12f, zwt_max_ - zwt_min_));
				}
				if (use_dcurv_) {
					row.push_back(
						(dcurv_[static_cast<std::size_t>(node_id)] - dcurv_min_) /
						std::max(1e-12f, dcurv_max_ - dcurv_min_));
				}

				const float y = response_[static_cast<std::size_t>(node_id)];
				const float w = scaled_weight_workspace_[static_cast<std::size_t>(position)];
				for (int j = 0; j < n_var; ++j) {
					XtY[static_cast<std::size_t>(j)] += w * row[static_cast<std::size_t>(j)] * y;
					for (int l = 0; l < n_var; ++l) {
						XtX[static_cast<std::size_t>(j)][static_cast<std::size_t>(l)] +=
							w * row[static_cast<std::size_t>(j)] * row[static_cast<std::size_t>(l)];
					}
				}
			}

			std::vector<std::vector<float>> inverse;
			if (!invert_matrix(XtX, inverse)) return false;

			beta.assign(static_cast<std::size_t>(n_var), 0.0f);
			for (int j = 0; j < n_var; ++j) {
				for (int l = 0; l < n_var; ++l) {
					beta[static_cast<std::size_t>(j)] +=
						inverse[static_cast<std::size_t>(j)][static_cast<std::size_t>(l)] *
						XtY[static_cast<std::size_t>(l)];
				}
			}
			return true;
		}

	private:
		/**
		 * @brief Initializes normalization metadata and local redundancy neighborhoods.
		 */
		void initialize()
		{
			observation_count_ = static_cast<int>(indices_.size());
			if (observation_count_ < 32 || (!use_zwt_ && !use_dcurv_)) return;

			const int loo_count = observation_count_ - 1;
			const int full_k = std::max(1, std::min(8, observation_count_ / 10));
			const int loo_k = std::max(1, std::min(8, loo_count / 10));
			if (full_k != loo_k) return;
			K_ = full_k;

			compute_extrema_metadata();
			base_raw_weight_.assign(static_cast<std::size_t>(observation_count_), 1.0f);
			affected_by_exclusion_.assign(
				static_cast<std::size_t>(observation_count_), {});
			override_raw_.assign(static_cast<std::size_t>(observation_count_), 1.0f);
			override_generation_.assign(static_cast<std::size_t>(observation_count_), 0u);
			scaled_weight_workspace_.assign(static_cast<std::size_t>(observation_count_), 1.0f);

			if (use_zwt_ && use_dcurv_) {
				initialize_two_dimensional();
			}
			else {
				initialize_one_dimensional();
			}

			if (!enabled_) return;
			raw_weight_order_.reserve(static_cast<std::size_t>(observation_count_));
			for (int position = 0; position < observation_count_; ++position) {
				raw_weight_order_.emplace_back(
					base_raw_weight_[static_cast<std::size_t>(position)], position);
			}
			std::sort(raw_weight_order_.begin(), raw_weight_order_.end());
		}

		/**
		 * @brief Builds the two-predictor nearest-neighbor cache.
		 */
		void initialize_two_dimensional()
		{
			local_count_ = std::min(observation_count_, 2 * K_ + 1);
			const int loo_local_count = std::min(observation_count_ - 1, 2 * K_ + 1);
			if (local_count_ != loo_local_count || local_count_ >= observation_count_) return;

			z01_.reserve(static_cast<std::size_t>(observation_count_));
			d01_.reserve(static_cast<std::size_t>(observation_count_));
			const float zrng = std::max(1e-12f, zwt_max_ - zwt_min_);
			const float drng = std::max(1e-12f, dcurv_max_ - dcurv_min_);
			for (int node_id : indices_) {
				z01_.push_back((zwt_[static_cast<std::size_t>(node_id)] - zwt_min_) / zrng);
				d01_.push_back((dcurv_[static_cast<std::size_t>(node_id)] - dcurv_min_) / drng);
			}

			const int stored_count = local_count_ + 1;
			nearest_two_dimensional_.assign(
				static_cast<std::size_t>(observation_count_),
				std::vector<std::pair<float, int>>(static_cast<std::size_t>(stored_count)));
			std::vector<std::pair<float, int>> local_dist(
				static_cast<std::size_t>(observation_count_));

			for (int k = 0; k < observation_count_; ++k) {
				for (int s = 0; s < observation_count_; ++s) {
					const float dz = z01_[static_cast<std::size_t>(s)] - z01_[static_cast<std::size_t>(k)];
					const float dd = d01_[static_cast<std::size_t>(s)] - d01_[static_cast<std::size_t>(k)];
					local_dist[static_cast<std::size_t>(s)] =
					{ std::sqrt(dz * dz + dd * dd), s };
				}
				std::partial_sort(
					local_dist.begin(),
					local_dist.begin() + stored_count,
					local_dist.end());

				float acc = 0.0f;
				float norm = 0.0f;
				for (int s = 0; s < stored_count; ++s) {
					nearest_two_dimensional_[static_cast<std::size_t>(k)][static_cast<std::size_t>(s)] =
						local_dist[static_cast<std::size_t>(s)];
					if (s < local_count_) {
						const float distance = local_dist[static_cast<std::size_t>(s)].first;
						const float kernel = std::max(0.0f, 1.0f - distance);
						acc += kernel;
						norm += 1.0f;
						affected_by_exclusion_[static_cast<std::size_t>(
							local_dist[static_cast<std::size_t>(s)].second)].push_back(k);
					}
				}
				const float density = (norm > 0.0f ? acc / norm : 1.0f);
				base_raw_weight_[static_cast<std::size_t>(k)] = 1.0f / (1.0f + density);
			}
			enabled_ = true;
		}

		/**
		 * @brief Builds the one-predictor sorted-window cache.
		 */
		void initialize_one_dimensional()
		{
			axis01_.resize(static_cast<std::size_t>(observation_count_));
			if (use_dcurv_) {
				const float range = std::max(1e-12f, dcurv_max_ - dcurv_min_);
				for (int position = 0; position < observation_count_; ++position) {
					const int node_id = indices_[static_cast<std::size_t>(position)];
					axis01_[static_cast<std::size_t>(position)] =
						(dcurv_[static_cast<std::size_t>(node_id)] - dcurv_min_) / range;
				}
			}
			else {
				const float range = std::max(1e-12f, zwt_max_ - zwt_min_);
				for (int position = 0; position < observation_count_; ++position) {
					const int node_id = indices_[static_cast<std::size_t>(position)];
					axis01_[static_cast<std::size_t>(position)] =
						(zwt_[static_cast<std::size_t>(node_id)] - zwt_min_) / range;
				}
			}

			one_dimensional_order_.reserve(static_cast<std::size_t>(observation_count_));
			for (int position = 0; position < observation_count_; ++position) {
				one_dimensional_order_.emplace_back(
					axis01_[static_cast<std::size_t>(position)], position);
			}
			std::sort(one_dimensional_order_.begin(), one_dimensional_order_.end());
			one_dimensional_rank_.assign(static_cast<std::size_t>(observation_count_), -1);
			for (int rank = 0; rank < observation_count_; ++rank) {
				one_dimensional_rank_[static_cast<std::size_t>(
					one_dimensional_order_[static_cast<std::size_t>(rank)].second)] = rank;
			}

			for (int rank = 0; rank < observation_count_; ++rank) {
				const int position = one_dimensional_order_[static_cast<std::size_t>(rank)].second;
				const int left = std::max(0, rank - K_);
				const int right = std::min(observation_count_ - 1, rank + K_);
				float acc = 0.0f;
				float norm = 0.0f;
				for (int s = left; s <= right; ++s) {
					const int neighbor_position =
						one_dimensional_order_[static_cast<std::size_t>(s)].second;
					const float distance = std::abs(
						one_dimensional_order_[static_cast<std::size_t>(s)].first -
						one_dimensional_order_[static_cast<std::size_t>(rank)].first);
					const float kernel = std::max(0.0f, 1.0f - distance);
					acc += kernel;
					norm += 1.0f;
					affected_by_exclusion_[static_cast<std::size_t>(neighbor_position)].push_back(position);
				}
				const float density = (norm > 0.0f ? acc / norm : 1.0f);
				base_raw_weight_[static_cast<std::size_t>(position)] = 1.0f / (1.0f + density);
			}
			enabled_ = true;
		}

		/**
		 * @brief Computes exact full-subset extrema and their multiplicities.
		 */
		void compute_extrema_metadata()
		{
			auto compute = [&](const std::vector<float>& values,
				float& minimum, float& maximum, int& minimum_count, int& maximum_count) {
				minimum = std::numeric_limits<float>::infinity();
				maximum = -std::numeric_limits<float>::infinity();
				for (int node_id : indices_) {
					const float value = values[static_cast<std::size_t>(node_id)];
					minimum = std::min(minimum, value);
					maximum = std::max(maximum, value);
				}
				minimum_count = 0;
				maximum_count = 0;
				for (int node_id : indices_) {
					const float value = values[static_cast<std::size_t>(node_id)];
					if (value == minimum) ++minimum_count;
					if (value == maximum) ++maximum_count;
				}
			};

			if (use_zwt_) {
				compute(zwt_, zwt_actual_min_, zwt_actual_max_, zwt_min_count_, zwt_max_count_);
			}
			if (use_dcurv_) {
				compute(dcurv_, dcurv_actual_min_, dcurv_actual_max_, dcurv_min_count_, dcurv_max_count_);
			}
		}

		/**
		 * @brief Tests whether removing a value changes a predictor normalization range.
		 */
		static bool exclusion_changes_range(
			const float value,
			const float minimum,
			const float maximum,
			const int minimum_count,
			const int maximum_count)
		{
			return (value == minimum && minimum_count == 1) ||
				(value == maximum && maximum_count == 1);
		}

		/**
		 * @brief Recomputes one raw redundancy weight after one observation is removed.
		 */
		float recompute_raw_weight(const int position, const int excluded_position) const
		{
			if (use_zwt_ && use_dcurv_) {
				float acc = 0.0f;
				float norm = 0.0f;
				int accepted = 0;
				for (const auto& item : nearest_two_dimensional_[static_cast<std::size_t>(position)]) {
					if (item.second == excluded_position) continue;
					const float kernel = std::max(0.0f, 1.0f - item.first);
					acc += kernel;
					norm += 1.0f;
					if (++accepted == local_count_) break;
				}
				const float density = (norm > 0.0f ? acc / norm : 1.0f);
				return 1.0f / (1.0f + density);
			}

			const int excluded_rank =
				one_dimensional_rank_[static_cast<std::size_t>(excluded_position)];
			const int full_rank = one_dimensional_rank_[static_cast<std::size_t>(position)];
			const int loo_rank = full_rank - (excluded_rank < full_rank ? 1 : 0);
			const int loo_count = observation_count_ - 1;
			const int left = std::max(0, loo_rank - K_);
			const int right = std::min(loo_count - 1, loo_rank + K_);

			float acc = 0.0f;
			float norm = 0.0f;
			const float target = axis01_[static_cast<std::size_t>(position)];
			for (int loo_sorted_position = left; loo_sorted_position <= right; ++loo_sorted_position) {
				const int full_sorted_position =
					loo_sorted_position >= excluded_rank
					? loo_sorted_position + 1
					: loo_sorted_position;
				const float distance = std::abs(
					one_dimensional_order_[static_cast<std::size_t>(full_sorted_position)].first - target);
				const float kernel = std::max(0.0f, 1.0f - distance);
				acc += kernel;
				norm += 1.0f;
			}
			const float density = (norm > 0.0f ? acc / norm : 1.0f);
			return 1.0f / (1.0f + density);
		}

		/**
		 * @brief Starts a new sparse raw-weight override generation.
		 */
		void begin_override_generation()
		{
			++override_generation_id_;
			if (override_generation_id_ == 0u) {
				std::fill(override_generation_.begin(), override_generation_.end(), 0u);
				override_generation_id_ = 1u;
			}
		}

		/**
		 * @brief Tests whether one observation has a raw-weight override.
		 */
		bool has_override(const int position) const
		{
			return override_generation_[static_cast<std::size_t>(position)] == override_generation_id_;
		}

		/**
		 * @brief Returns the active raw redundancy weight of one observation.
		 */
		float raw_weight_at(const int position) const
		{
			return has_override(position)
				? override_raw_[static_cast<std::size_t>(position)]
				: base_raw_weight_[static_cast<std::size_t>(position)];
		}

		const std::vector<int>& indices_;
		const std::vector<float>& zwt_;
		const std::vector<float>& dcurv_;
		const std::vector<float>& response_;
		const std::vector<std::uint8_t>& is_spring_node_;
		bool use_zwt_ = false;
		bool use_dcurv_ = false;
		float zwt_min_ = 0.0f;
		float zwt_max_ = 1.0f;
		float dcurv_min_ = 0.0f;
		float dcurv_max_ = 1.0f;
		int observation_count_ = 0;
		int K_ = 0;
		int local_count_ = 0;
		bool enabled_ = false;

		float zwt_actual_min_ = 0.0f;
		float zwt_actual_max_ = 0.0f;
		float dcurv_actual_min_ = 0.0f;
		float dcurv_actual_max_ = 0.0f;
		int zwt_min_count_ = 0;
		int zwt_max_count_ = 0;
		int dcurv_min_count_ = 0;
		int dcurv_max_count_ = 0;

		std::vector<float> z01_;
		std::vector<float> d01_;
		std::vector<float> axis01_;
		std::vector<std::vector<std::pair<float, int>>> nearest_two_dimensional_;
		std::vector<std::pair<float, int>> one_dimensional_order_;
		std::vector<int> one_dimensional_rank_;
		std::vector<std::vector<int>> affected_by_exclusion_;
		std::vector<float> base_raw_weight_;
		std::vector<std::pair<float, int>> raw_weight_order_;
		std::vector<float> override_raw_;
		std::vector<std::uint32_t> override_generation_;
		std::uint32_t override_generation_id_ = 0u;
		std::vector<float> scaled_weight_workspace_;
	};
}


namespace {
	/**
	 * @brief Returns the two closest entries of a uniquely valued sorted distribution.
	 *
	 * The ordering criterion is identical to `find_closest_values`: absolute
	 * distance to the query first, then the smaller property value on an exact
	 * distance tie. Only the two bracketing entries and their immediate neighbors
	 * can be among the two closest values on a sorted one-dimensional support.
	 *
	 * @param sorted_values Pairs of (value, original index), sorted by value.
	 * @param query Query value.
	 * @return Two closest (value, original index) pairs in the legacy ordering.
	 */
	static std::pair<std::pair<float, std::size_t>, std::pair<float, std::size_t>>
		find_two_closest_unique_sorted_values(
			const std::vector<std::pair<float, std::size_t>>& sorted_values,
			const float query)
	{
		auto lower = std::lower_bound(
			sorted_values.begin(), sorted_values.end(), query,
			[](const std::pair<float, std::size_t>& item, const float value) {
			return item.first < value;
		});
		const int center = static_cast<int>(lower - sorted_values.begin());
		const int first = std::max(0, center - 2);
		const int last = std::min(static_cast<int>(sorted_values.size()), center + 2);

		std::array<std::pair<float, std::size_t>, 4> candidates{};
		int candidate_count = 0;
		for (int index = first; index < last; ++index) {
			candidates[static_cast<std::size_t>(candidate_count++)] =
				sorted_values[static_cast<std::size_t>(index)];
		}
		std::sort(
			candidates.begin(), candidates.begin() + candidate_count,
			[query](const auto& lhs, const auto& rhs) {
			const float lhs_diff = std::abs(lhs.first - query);
			const float rhs_diff = std::abs(rhs.first - query);
			if (lhs_diff == rhs_diff) return lhs.first < rhs.first;
			return lhs_diff < rhs_diff;
		});
		return { candidates[0], candidates[1] };
	}

	/**
	 * @brief Performs the conditioning-data normal-score transform without repeated full sorts.
	 *
	 * The legacy implementation calls `find_closest_values` independently for each
	 * conditioning datum, which rebuilds and sorts the complete simulation
	 * distribution every time. When distribution values are unique, the same two
	 * values are obtained by one global value sort followed by logarithmic searches.
	 * If exact duplicate values are present, this helper deliberately falls back to
	 * the legacy routine because `std::sort` does not define the relative order of
	 * comparator-equivalent duplicates and the legacy index permutation is observable
	 * in its current quantile lookup.
	 *
	 * @param data_vector Conditioning values, including optional -99999 NDV entries.
	 * @param discrete_distribution Simulation marginal distribution.
	 * @param mean Target Gaussian mean.
	 * @param stddev Target Gaussian standard deviation.
	 * @return Normal-score transformed conditioning vector.
	 */
	static std::vector<float> nst_data_with_nodata_fast_equivalent(
		const std::vector<float>& data_vector,
		const std::vector<float>& discrete_distribution,
		const float mean,
		const float stddev)
	{
		std::vector<float> filtered_distribution;
		filtered_distribution.reserve(discrete_distribution.size());
		for (const float value : discrete_distribution) {
			if (std::abs(value - (-99999.0f)) > 1e-12f) {
				if (!std::isfinite(value)) {
					return nst_data_with_nodata(data_vector, discrete_distribution, mean, stddev);
				}
				filtered_distribution.push_back(value);
			}
		}
		if (filtered_distribution.size() < 2u) {
			return nst_data_with_nodata(data_vector, discrete_distribution, mean, stddev);
		}
		for (const float value : data_vector) {
			if (std::abs(value - (-99999.0f)) > 1e-12f && !std::isfinite(value)) {
				return nst_data_with_nodata(data_vector, discrete_distribution, mean, stddev);
			}
		}

		std::vector<std::pair<float, std::size_t>> sorted_values;
		sorted_values.reserve(filtered_distribution.size());
		for (std::size_t index = 0; index < filtered_distribution.size(); ++index) {
			sorted_values.emplace_back(filtered_distribution[index], index);
		}
		std::sort(
			sorted_values.begin(), sorted_values.end(),
			[](const auto& lhs, const auto& rhs) { return lhs.first < rhs.first; });

		for (std::size_t index = 1; index < sorted_values.size(); ++index) {
			if (sorted_values[index].first == sorted_values[index - 1].first) {
				return nst_data_with_nodata(data_vector, discrete_distribution, mean, stddev);
			}
		}

		std::vector<std::size_t> ranks(filtered_distribution.size());
		std::iota(ranks.begin(), ranks.end(), 0);
		std::sort(
			ranks.begin(), ranks.end(),
			[&](const std::size_t i, const std::size_t j) {
			return filtered_distribution[i] > filtered_distribution[j];
		});

		std::vector<float> quantiles(filtered_distribution.size());
		for (std::size_t rank = 0; rank < filtered_distribution.size(); ++rank) {
			quantiles[ranks[rank]] =
				(static_cast<float>(rank) + 0.5f) /
				static_cast<float>(filtered_distribution.size());
		}

		std::vector<float> transformed_values(data_vector.size());
		for (std::size_t index = 0; index < data_vector.size(); ++index) {
			const float data = data_vector[index];
			if (std::abs(data - (-99999.0f)) <= 1e-12f) {
				transformed_values[index] = -99999.0f;
				continue;
			}

			const auto closest = find_two_closest_unique_sorted_values(sorted_values, data);
			const float q1 = quantiles[ranks[closest.first.second]];
			const float q2 = quantiles[ranks[closest.second.second]];
			float q = 0.0f;
			const bool has_two_neighbors =
				std::abs(closest.first.first - data) >= 1e-4f &&
				!((closest.first.first > data && closest.second.first > data) ||
					(closest.first.first < data && closest.second.first < data));

			if (has_two_neighbors) {
				q = interpolate(
					q1, closest.first.first, q2, closest.second.first, data);
			}
			else if (closest.first.first < data) {
				q = interpolate(
					q1, closest.first.first, q2, closest.second.first, data);
			}
			else if (closest.first.first > data) {
				q = interpolate(
					q1, closest.first.first, q2, closest.second.first, data);
			}
			else {
				q = q1;
			}
			transformed_values[index] = mean + stddev * inverse_normal_cdf(1.0f - q);
		}
		return transformed_values;
	}

	/**
	 * @brief Back-transforms Gaussian values using one ranked marginal-distribution table.
	 *
	 * Quantiles of the initial empirical distribution are uniformly spaced rank
	 * midpoints. The legacy routine nevertheless sorts the complete quantile table
	 * for every simulated node to recover the two closest ranks. This implementation
	 * constructs the same ranked property values and midpoint quantiles once, then
	 * searches only the local bracketing ranks. The interpolation operands and their
	 * order are identical to the legacy implementation, including extrapolation.
	 *
	 * @param gaussian_distribution Gaussian-space simulated values.
	 * @param discrete_distribution Initial property-space marginal distribution.
	 * @return Property-space back-transformed values.
	 */
	static std::vector<float> back_transform_fast_equivalent(
		const std::vector<float>& gaussian_distribution,
		const std::vector<float>& discrete_distribution)
	{
		std::vector<float> filtered_distribution;
		filtered_distribution.reserve(discrete_distribution.size());
		for (const float value : discrete_distribution) {
			if (std::abs(value - (-99999.0f)) > 1e-12f) {
				if (!std::isfinite(value)) {
					return back_transform(gaussian_distribution, discrete_distribution);
				}
				filtered_distribution.push_back(value);
			}
		}
		if (filtered_distribution.size() < 2u) {
			return back_transform(gaussian_distribution, discrete_distribution);
		}

		std::vector<std::size_t> ranks(filtered_distribution.size());
		std::iota(ranks.begin(), ranks.end(), 0);
		std::sort(
			ranks.begin(), ranks.end(),
			[&](const std::size_t i, const std::size_t j) {
			return filtered_distribution[i] < filtered_distribution[j];
		});

		std::vector<float> sorted_property(filtered_distribution.size());
		std::vector<float> rank_quantiles(filtered_distribution.size());
		for (std::size_t rank = 0; rank < filtered_distribution.size(); ++rank) {
			sorted_property[rank] = filtered_distribution[ranks[rank]];
			rank_quantiles[rank] =
				(static_cast<float>(rank) + 0.5f) /
				static_cast<float>(filtered_distribution.size());
		}

		std::vector<float> back_transformed(gaussian_distribution.size());
		for (std::size_t index = 0; index < gaussian_distribution.size(); ++index) {
			const float gaussian_value = gaussian_distribution[index];
			if (std::abs(gaussian_value - (-99999.0f)) <= 1e-12f) {
				back_transformed[index] = -99999.0f;
				continue;
			}

			const float tested_quantile = normal_cdf(gaussian_value);
			auto lower = std::lower_bound(
				rank_quantiles.begin(), rank_quantiles.end(), tested_quantile);
			const int center = static_cast<int>(lower - rank_quantiles.begin());
			const int first = std::max(0, center - 2);
			const int last = std::min(static_cast<int>(rank_quantiles.size()), center + 2);

			std::array<std::pair<float, int>, 4> candidates{};
			int candidate_count = 0;
			for (int rank = first; rank < last; ++rank) {
				candidates[static_cast<std::size_t>(candidate_count++)] =
				{ rank_quantiles[static_cast<std::size_t>(rank)], rank };
			}
			std::sort(
				candidates.begin(), candidates.begin() + candidate_count,
				[tested_quantile](const auto& lhs, const auto& rhs) {
				const float lhs_diff = std::abs(lhs.first - tested_quantile);
				const float rhs_diff = std::abs(rhs.first - tested_quantile);
				if (lhs_diff == rhs_diff) return lhs.first < rhs.first;
				return lhs_diff < rhs_diff;
			});

			const float q1 = candidates[0].first;
			const float q2 = candidates[1].first;
			const float value1 = sorted_property[static_cast<std::size_t>(candidates[0].second)];
			const float value2 = sorted_property[static_cast<std::size_t>(candidates[1].second)];
			const bool has_two_neighbors =
				std::abs(q1 - tested_quantile) >= 1e-4f &&
				!((q1 > tested_quantile && q2 > tested_quantile) ||
					(q1 < tested_quantile && q2 < tested_quantile));

			if (has_two_neighbors) {
				back_transformed[index] = interpolate(
					value1, q1, value2, q2, tested_quantile);
			}
			else if (q1 < tested_quantile) {
				back_transformed[index] = interpolate(
					value1, q1, value2, q2, tested_quantile);
			}
			else if (q1 > tested_quantile) {
				back_transformed[index] = interpolate(
					value1, q1, value2, q2, tested_quantile);
			}
			else {
				back_transformed[index] = value1;
			}
		}
		return back_transformed;
	}
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

	auto& workspace = GEOSTATS_DIJKSTRA_WORKSPACE;
	workspace.begin(static_cast<std::size_t>(N));

	struct Q { float d; int i; bool operator>(const Q& o) const { return d > o.d; } };
	std::priority_queue<Q, std::vector<Q>, std::greater<Q>> pq;

	auto edge_len = [&](int u, int v)->float {
		const Vector3& pu = curve->nodes[u].p;
		const Vector3& pv = curve->nodes[v].p;
		return KarstNSim::magnitude(pu - pv);
	};

	workspace.set_distance(src, 0.f);
	pq.push({ 0.f, src });

	// Preserve the original target counting semantics, including invalid or
	// duplicate target entries, while avoiding an O(N) target-array reset.
	for (int t : targets) {
		if (t >= 0 && t < N) workspace.mark_target(t);
	}
	int remaining = (int)targets.size();

	while (!pq.empty()) {
		auto [du, u] = pq.top(); pq.pop();
		if (workspace.is_settled(u)) continue;
		workspace.mark_settled(u);

		if (du > range_cap) break;               // all further keys >= du
		if (workspace.is_target(u)) {            // settled a target
			if (--remaining == 0) break;          // all found
		}

		for (const auto& c : curve->nodes[u].connections) {
			const int v = c.destindex;
			if (v < 0 || v >= N) continue;
			if (workspace.is_settled(v)) continue;
			const float w = edge_len(u, v);
			if (!std::isfinite(w) || w < 0.f) continue;
			const float alt = du + w;
			if (alt >= workspace.get_distance(v) || alt > range_cap) continue;
			workspace.set_distance(v, alt);
			pq.push({ alt, v });
		}
	}

	std::vector<float> out; out.reserve(targets.size());
	for (int t : targets) out.push_back((t >= 0 && t < N) ? workspace.get_distance(t) : INF);
	return out;
}

/**
 * @brief Performs simple kriging using graph distances computed on the fly.
 *
 * The local covariance system is assembled from shortest-path distances and
 * solved directly by a Cholesky decomposition in double precision. No global
 * distance matrix is constructed. The known mean is the empirical mean of the
 * Gaussian kriging distribution, consistently with the former simple-kriging
 * implementation.
 *
 * When the neighborhood is empty, or when the local system cannot be solved,
 * the function returns one unconditional draw from the Gaussian distribution.
 */
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
	// If the neighborhood is empty, perform one unconditional empirical draw.
	if (neighborhood.empty()) {
		val_estimation = select_random_element(*kriging_distribution);
		var_estimation = static_cast<float>(active_zero_lag_variance(vario_sill));
		return;
	}

	auto unconditional_fallback = [&](const std::string& reason) {
		const Vector3& point = curve->nodes.at(current_node_index).p;
		std::cerr
			<< "[geostats][simple_kriging] Unconditional fallback at skeleton node "
			<< current_node_index << " (X=" << point.x
			<< ", Y=" << point.y << ", Z=" << point.z << "): "
			<< reason << "." << std::endl;
		val_estimation = select_random_element(*kriging_distribution);
		var_estimation = std::numeric_limits<float>::quiet_NaN();
	};

	const int K = (int)neighborhood.size();
	// 1) Distances current -> neighbors (shortest path, truncated at range_cap)
	std::vector<int> targets = neighborhood;
	std::vector<float> d_cur_to_nb = dijkstra_to_targets_truncated(
		current_node_index, curve, targets, range_cap
	);

	// 2) Pairwise neighbor-neighbor distances. Two nodes that are each at most R
	//    from the simulated node can be separated by up to 2R through that node.
	const float pairwise_range_cap =
		(range_cap <= std::numeric_limits<float>::max() / 2.0f)
		? 2.0f * range_cap
		: std::numeric_limits<float>::max();
	std::vector<std::vector<float>> d_nb_to_nb(K, std::vector<float>(K, 0.f));
	for (int i = 0; i < K; ++i) {
		std::vector<float> di = dijkstra_to_targets_truncated(
			neighborhood[i], curve, targets, pairwise_range_cap);
		for (int j = 0; j < K; ++j) d_nb_to_nb[i][j] = di[j];
	}

	// 3) Build the simple-kriging covariance system C * lambda = c0.
	const double zero_lag_variance = active_zero_lag_variance(vario_sill);
	auto covariance = [&](float distance)->double {
		if (distance <= 0.0f) return zero_lag_variance;
		const float semivariance = variogram_value(
			distance, vario_sill, vario_nugget, vario_range, vario_model);
		return zero_lag_variance - static_cast<double>(semivariance);
	};

	std::vector<std::vector<double>> covariance_matrix(
		K, std::vector<double>(K, 0.0));
	std::vector<double> covariance_to_target(K, 0.0);

	for (int i = 0; i < K; ++i) {
		for (int j = 0; j < K; ++j) {
			const float distance = d_nb_to_nb[i][j];
			// An unreachable pair is treated as uncorrelated, consistently with
			// the previous use of the sill as its semivariance.
			covariance_matrix[i][j] =
				std::isfinite(distance) ? covariance(distance) : 0.0;
		}
		const float distance_to_target = d_cur_to_nb[i];
		covariance_to_target[i] =
			std::isfinite(distance_to_target) ? covariance(distance_to_target) : 0.0;
	}

	// 4) Solve directly; do not form the inverse covariance matrix.
	std::vector<double> weights;
	std::string failure_reason;
	if (!solve_spd_cholesky(
		covariance_matrix, covariance_to_target, weights, failure_reason)) {
		unconditional_fallback(failure_reason);
		return;
	}

	// 5) Simple-kriging estimate with known mean m:
	//    Z*(u0) = m + sum_i lambda_i [Z(ui) - m].
	const double known_mean = std::accumulate(
		kriging_distribution->begin(), kriging_distribution->end(), 0.0) /
		static_cast<double>(kriging_distribution->size());
	double estimate = known_mean;
	for (int i = 0; i < K; ++i) {
		estimate += weights[i] *
			(static_cast<double>(node_values[neighborhood[i]]) - known_mean);
	}
	if (!std::isfinite(estimate)) {
		unconditional_fallback("kriging estimate is non-finite");
		return;
	}

	// 6) Simple-kriging variance: sigma_K^2 = C(0) - lambda^T c0.
	double variance = zero_lag_variance;
	for (int i = 0; i < K; ++i) {
		variance -= weights[i] * covariance_to_target[i];
	}
	if (!std::isfinite(variance)) {
		unconditional_fallback("kriging variance is non-finite");
		return;
	}

	val_estimation = static_cast<float>(estimate);
	var_estimation = static_cast<float>(std::max(0.0, variance));
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
		const float supplied_sill = vario_sill;
		float C = static_cast<float>(active_zero_lag_variance(vario_sill));
		float C0 = vario_nugget;

		// Fill matrix of variogram values
		for (int j = 0; j < nb_neigh + 1; ++j) {
			float h = local_mat_distance[0][j];
			if (h == 0.) { // It should not be the case
				local_mat_vario_values[0][j] = 0.;
			}
			else {
				local_mat_vario_values[0][j] = variogram_value(h, supplied_sill, C0, a, vario_model);
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
					local_mat_vario_values[i][j] = variogram_value(h, supplied_sill, C0, a, vario_model);
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
		// Linear scan is deliberately retained to preserve the exact tie-breaking
		// behavior of the previous implementation.
		int best = -1;
		float best2 = std::numeric_limits<float>::infinity();
		for (int i = 0; i < (int)sk->nodes.size(); ++i) {
			const auto& Q = sk->nodes[i].p;
			float d2 = sqr(P.x - Q.x) + sqr(P.y - Q.y) + sqr(P.z - Q.z);
			if (d2 < best2) { best2 = d2; best = i; }
		}
		return best;
	}

	/**
	 * @brief Stores the nearest skeleton node and elevation of one input spring.
	 */
	struct MappedSpring {
		int node_id;
		float spring_z;
	};

	/**
	 * @brief Maps every input spring to its nearest skeleton node once.
	 *
	 * The former implementation repeated the same O(N) nearest-node scan for every
	 * spring in every connected component. The mapping is independent of the
	 * component being processed, so computing it once preserves exactly the same
	 * node selection while removing redundant work.
	 *
	 * @param sk Skeleton graph.
	 * @param springs_xyz Input spring coordinates.
	 * @return Spring-to-node mapping in the original spring order.
	 */
	static std::vector<MappedSpring> map_springs_to_nodes(
		const KarstNSim::KarsticSkeleton* sk,
		const std::vector<Vector3>& springs_xyz)
	{
		std::vector<MappedSpring> mapped;
		mapped.reserve(springs_xyz.size());
		for (const Vector3& spring : springs_xyz) {
			mapped.push_back({ nearest_node_index(sk, spring), spring.z });
		}
		return mapped;
	}

	/**
	 * @brief Selects source nodes for one connected component from pre-mapped springs.
	 *
	 * @param comp Connected-component label of every skeleton node.
	 * @param cid Component being processed.
	 * @param mapped_springs Precomputed nearest-node mapping of all input springs.
	 * @param z_window Elevation window above the lowest spring.
	 * @return Deduplicated source-node indices in ascending node-index order.
	 */
	static std::vector<int> pick_component_sources_from_mapped_springs(
		const std::vector<int>& comp,
		int cid,
		const std::vector<MappedSpring>& mapped_springs,
		float z_window = 40.f)
	{
		std::vector<MappedSpring> candidates;
		candidates.reserve(mapped_springs.size());

		for (const MappedSpring& spring : mapped_springs) {
			if (spring.node_id >= 0 && comp[spring.node_id] == cid) {
				candidates.push_back(spring);
			}
		}
		if (candidates.empty()) return {};

		float zmin = candidates[0].spring_z;
		for (const MappedSpring& candidate : candidates) {
			zmin = std::min(zmin, candidate.spring_z);
		}

		std::vector<int> sources_nodes;
		sources_nodes.reserve(candidates.size());
		for (const MappedSpring& candidate : candidates) {
			if (candidate.spring_z <= zmin + z_window) {
				sources_nodes.push_back(candidate.node_id);
			}
		}

		std::sort(sources_nodes.begin(), sources_nodes.end());
		sources_nodes.erase(
			std::unique(sources_nodes.begin(), sources_nodes.end()),
			sources_nodes.end());
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

	/**
	 * @brief Upstream DAG edge with a compact union identifier and cached length.
	 */
	struct UpstreamParentEdge {
		int upstream_node;
		int edge_id;
		float length;
	};

	static int build_upstream_parents(
		const KarstNSim::KarsticSkeleton* sk,
		const std::vector<int>& comp, int cid,
		const std::vector<int>& component_nodes,
		const std::vector<float>& dist,
		std::vector<std::vector<UpstreamParentEdge>>& parents)
	{
		const int N = (int)sk->nodes.size();
		parents.assign(N, {});

		// Assign a compact identifier to each undirected DAG edge once. Duplicate
		// adjacency entries, if any, deliberately receive the same identifier so
		// that the former unordered-set edge-union semantics are preserved.
		std::unordered_map<uint64_t, int> edge_ids;
		edge_ids.reserve(std::max<std::size_t>(1024, component_nodes.size() * 2));
		int next_edge_id = 0;

		for (int v : component_nodes) {
			for (const auto& c : sk->nodes[v].connections) {
				const int u = c.destindex;
				if (u < 0 || u >= N) continue;
				if (comp[u] != cid) continue;
				if (dist[u] > dist[v] + DCURV_EPS) {
					const uint64_t key = make_edge_key(u, v);
					auto insertion = edge_ids.emplace(key, next_edge_id);
					if (insertion.second) ++next_edge_id;
					parents[v].push_back({
						u,
						insertion.first->second,
						edge_length(sk, u, v)
						});
				}
			}
		}
		return next_edge_id;
	}

	// --- Exact upstream union per target set ---------------------------------
	// Each edge is counted at most once per target by using an unordered_set of undirected edge keys.

	static std::vector<float> exact_dcurv_for_targets(
		const KarstNSim::KarsticSkeleton* curve,
		const std::vector<int>& comp, int cid,
		const std::vector<std::vector<UpstreamParentEdge>>& parents,
		const int edge_count,
		const std::vector<int>& targets)
	{
		std::vector<float> out;
		out.reserve(targets.size());

		std::vector<int> stack;
		stack.reserve(1024);

		// Reusable generation-stamped arrays replace one unordered_set allocation
		// and thousands of hash operations per target. The DFS stack order and the
		// order in which edge lengths are added are unchanged, so floating-point
		// accumulation follows the same sequence as in the previous implementation.
		std::vector<std::uint32_t> seen_edge_generation(
			static_cast<std::size_t>(std::max(0, edge_count)), 0u);
		std::vector<std::uint32_t> expanded_node_generation(
			curve->nodes.size(), 0u);
		std::uint32_t generation = 0u;

		for (int t : targets) {
			++generation;
			if (generation == 0u) {
				std::fill(seen_edge_generation.begin(), seen_edge_generation.end(), 0u);
				std::fill(expanded_node_generation.begin(), expanded_node_generation.end(), 0u);
				generation = 1u;
			}

			if (t < 0 || t >= (int)curve->nodes.size() || comp[t] != cid) {
				out.push_back(0.0f);
				continue;
			}

			float total = 0.0f;
			stack.clear();
			stack.push_back(t);

			// Non-recursive DFS up the parents DAG. A node reached more than once can
			// be skipped after its first expansion because all of its parent edges were
			// already examined atomically during that expansion.
			while (!stack.empty()) {
				const int v = stack.back();
				stack.pop_back();

				if (expanded_node_generation[static_cast<std::size_t>(v)] == generation) {
					continue;
				}
				expanded_node_generation[static_cast<std::size_t>(v)] = generation;

				for (const UpstreamParentEdge& parent : parents[v]) {
					const std::size_t edge_index = static_cast<std::size_t>(parent.edge_id);
					if (seen_edge_generation[edge_index] != generation) {
						seen_edge_generation[edge_index] = generation;
						total += parent.length;
						stack.push_back(parent.upstream_node);
					}
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
	int ncomp = 0;
	for (int x : comp) ncomp = std::max(ncomp, x + 1);

	// Cache component-node lists in increasing global node order. Besides avoiding
	// repeated O(N) scans for every component, this preserves the exact iteration
	// order used by the former implementation.
	std::vector<std::vector<int>> component_nodes(static_cast<std::size_t>(ncomp));
	for (int v = 0; v < N; ++v) {
		if (comp[v] >= 0) component_nodes[static_cast<std::size_t>(comp[v])].push_back(v);
	}

	// The nearest skeleton node of a spring is independent of the component loop.
	// Compute the identical linear-scan mapping once instead of once per component.
	const std::vector<MappedSpring> mapped_springs =
		map_springs_to_nodes(curve, springs_xyz);

	// 2) Process each component independently
	for (int cid = 0; cid < ncomp; ++cid) {
		const std::vector<int>& nodes_in_component =
			component_nodes[static_cast<std::size_t>(cid)];

		// 2.1) Select sources among springs located in this component:
		//      lowest-Z spring + any other spring within +40 m in Z.
		std::vector<int> sources = pick_component_sources_from_mapped_springs(
			comp, cid, mapped_springs, 40.f);

		if (sources.empty()) {
			continue;
		}

		// 2.2) Multi-source Dijkstra distances
		std::vector<float> dist = dijkstra_multi_source(curve, comp, cid, sources, nullptr);

		// 2.3) Upstream DAG (multi-parents). Parent edges receive compact IDs and
		//      cache their geometric length once for all target traversals.
		std::vector<std::vector<UpstreamParentEdge>> parents;
		const int upstream_edge_count = build_upstream_parents(
			curve, comp, cid, nodes_in_component, dist, parents);

		// Compute total undirected edge length of this component using the same
		// node and connection order as before.
		double comp_total_len = 0.0;
		for (int v : nodes_in_component) {
			for (const auto& c : curve->nodes[v].connections) {
				const int u = c.destindex;
				if (u < 0 || u >= N || comp[u] != cid) continue;
				if (u < v) comp_total_len += edge_length(curve, u, v);
			}
		}

		// 2.4) Exact upstream edge-union length for all nodes in the component.
		// `nodes_in_component` has the same increasing order as the former target
		// construction and can therefore be used directly.
		std::vector<float> dcurv_comp = exact_dcurv_for_targets(
			curve, comp, cid, parents, upstream_edge_count, nodes_in_component);

		for (std::size_t k = 0; k < nodes_in_component.size(); ++k) {
			const int v = nodes_in_component[k];
			dcurv[v] = dcurv_comp[k];
			if (dcurv[v] > comp_total_len + 1e-3f) {
				dcurv[v] = (float)comp_total_len;
			}
		}
	}

	// 3) Stats retained for strict behavioral parity with the previous code.
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
//     * local-density weighting in the joint normalized predictor space when
//       zwt and dcurv are both active, or along the single active predictor
// - Class-balance weighting:
//     * Within the current fitting subset, split total weight mass 50/50
//       between samples located at spring positions (SPRINGS) and all others (OTHERS).
//     * Uses exact proximity test with tolerance EPS = 1e-2 on (x,y,z),
//       via KarstNSim::magnitude(P - S) <= EPS.
// - Outlier rejection (hard 0/1):
//     * compute leave-one-out residuals whenever removal preserves the inlet
//       and outlet boundary classes and leaves a solvable regression
//     * robust scale via MAD (sigma = 1.4826 * median(|r_i|))
//     * reject testable observations if |r_i| / sigma > C_CUTOFF
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
	const std::vector<ConditioningDataRole>& conditioning_roles,
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

	auto conditioning_role_at = [&](const int node_id) -> ConditioningDataRole {
		if (node_id >= 0 && node_id < static_cast<int>(conditioning_roles.size())) {
			return conditioning_roles[node_id];
		}
		return ConditioningDataRole::None;
	};

	// ---- 1) Collect observed nodes (non-NDV eq_radius) ------------------------
	std::vector<int> valid_indices;
	valid_indices.reserve(N);

	int n_hard_data = 0;
	int n_excluded_by_radius_cap = 0;
	int n_candidate_inlets = 0;
	int n_candidate_outlets = 0;
	int n_candidate_waypoints = 0;

	for (int i = 0; i < N; ++i) {
		if (i >= static_cast<int>(eq_radius_values.size())) {
			continue;
		}

		const float r = eq_radius_values[i];
		const bool has_data =
			(std::abs(r - (-99999.0f)) > 1e-12f);

		if (!has_data) {
			continue;
		}

		++n_hard_data;

		const bool pass_radius_cap =
			(!K_EXT_DRIFT_ENABLE_RADIUS_CAP ||
				r <= K_RADIUS_MAX_FOR_REGRESSION);

		if (!pass_radius_cap) {
			++n_excluded_by_radius_cap;
			continue;
		}

		valid_indices.push_back(i);

		switch (conditioning_role_at(i)) {
		case ConditioningDataRole::Inlet:
			++n_candidate_inlets;
			break;
		case ConditioningDataRole::Outlet:
			++n_candidate_outlets;
			break;
		case ConditioningDataRole::Waypoint:
			++n_candidate_waypoints;
			break;
		default:
			break;
		}
	}

	if (K_EXT_DRIFT_ENABLE_DIAGNOSTIC_LOGS) {
		std::cout
			<< "[geostats][external_drift][diag] Hard conditioning data: "
			<< n_hard_data
			<< "; regression candidates: " << valid_indices.size()
			<< " (inlets=" << n_candidate_inlets
			<< ", outlets=" << n_candidate_outlets
			<< ", waypoints=" << n_candidate_waypoints << ")"
			<< "; excluded by radius cap: " << n_excluded_by_radius_cap;

		if (K_EXT_DRIFT_ENABLE_RADIUS_CAP) {
			std::cout << " (cap=" << K_RADIUS_MAX_FOR_REGRESSION << ")";
		}
		else {
			std::cout << " (cap disabled)";
		}

		std::cout << std::endl;
	}

	const int n_obs = static_cast<int>(valid_indices.size());
	if (n_obs == 0 || (!use_drift_zwt && !use_drift_curv)) {
		weights_out.assign(N, 0.0f);
		return drift;
	}

	// Precompute spring-class membership once. The previous implementation repeated
	// the same coordinate-to-spring proximity test inside every full and leave-one-
	// out regression fit, although the classification is subset-independent.
	const float SPR_EPS = 1e-2f;
	std::vector<std::uint8_t> is_spring_observation(static_cast<std::size_t>(N), 0u);
	for (int id : valid_indices) {
		const Vector3& point = curve->nodes[id].p;
		for (const Vector3& spring : springs_xyz) {
			if (KarstNSim::magnitude(point - spring) <= SPR_EPS) {
				is_spring_observation[static_cast<std::size_t>(id)] = 1u;
				break;
			}
		}
	}

	if (K_EXT_DRIFT_ENABLE_DIAGNOSTIC_LOGS) {
		int n_spring_class = 0;
		for (int id : valid_indices) {
			if (is_spring_observation[static_cast<std::size_t>(id)] != 0u) {
				++n_spring_class;
			}
		}

		std::cout
			<< "[geostats][external_drift][diag] Regression spring class: "
			<< n_spring_class
			<< " observations recognized geometrically, versus "
			<< n_candidate_outlets
			<< " observations tagged as outlets."
			<< std::endl;
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

	auto redundancy_weights = [&](const std::vector<int>& idxs,
		bool use_zwt, bool use_dcurv,
		const float& zwt_min, const float& zwt_max,
		const float& dcurv_min, const float& dcurv_max) {
		const int M = static_cast<int>(idxs.size());
		std::vector<float> w(M, 1.0f);
		if (M <= 1 || (!use_zwt && !use_dcurv)) return w;

		std::vector<float> z01;
		std::vector<float> d01;
		z01.reserve(M);
		d01.reserve(M);

		if (use_zwt) {
			const float rng = std::max(1e-12f, zwt_max - zwt_min);
			for (int id : idxs) z01.push_back((zwt[id] - zwt_min) / rng);
		}
		if (use_dcurv) {
			const float rng = std::max(1e-12f, dcurv_max - dcurv_min);
			for (int id : idxs) d01.push_back((dcurv[id] - dcurv_min) / rng);
		}

		const int K = std::max(1, std::min(8, M / 10));

		if (use_zwt && use_dcurv) {
			const int local_count = std::min(M, 2 * K + 1);
			std::vector<std::pair<float, int>> local_dist(static_cast<std::size_t>(M));

			for (int k = 0; k < M; ++k) {
				for (int s = 0; s < M; ++s) {
					const float dz = z01[s] - z01[k];
					const float dd = d01[s] - d01[k];
					local_dist[static_cast<std::size_t>(s)] =
					{ std::sqrt(dz * dz + dd * dd), s };
				}

				// Only the nearest `local_count` entries contribute to the density.
				// partial_sort returns exactly the same ordered prefix as a complete
				// lexicographic sort of (distance, index), while reducing the sorting
				// cost from O(M log M) to O(M log local_count), with local_count <= 17.
				if (local_count < M) {
					std::partial_sort(
						local_dist.begin(),
						local_dist.begin() + local_count,
						local_dist.end());
				}
				else {
					std::sort(local_dist.begin(), local_dist.end());
				}

				float acc = 0.0f, norm = 0.0f;
				for (int s = 0; s < local_count; ++s) {
					const float distance = local_dist[static_cast<std::size_t>(s)].first;
					const float kernel = std::max(0.0f, 1.0f - distance);
					acc += kernel;
					norm += 1.0f;
				}
				const float density = (norm > 0.0f ? acc / norm : 1.0f);
				w[k] = 1.0f / (1.0f + density);
			}
		}
		else {
			const std::vector<float>& axis01 = use_dcurv ? d01 : z01;
			std::vector<std::pair<float, int>> order;
			order.reserve(M);
			for (int k = 0; k < M; ++k) order.emplace_back(axis01[k], k);
			std::sort(order.begin(), order.end());

			for (int t = 0; t < M; ++t) {
				const int k0 = order[t].second;
				const int kL = std::max(0, t - K);
				const int kR = std::min(M - 1, t + K);
				float acc = 0.0f, norm = 0.0f;
				for (int s = kL; s <= kR; ++s) {
					const float distance = std::abs(order[s].first - order[t].first);
					const float kernel = std::max(0.0f, 1.0f - distance);
					acc += kernel;
					norm += 1.0f;
				}
				const float density = (norm > 0.0f ? acc / norm : 1.0f);
				w[k0] = 1.0f / (1.0f + density);
			}
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

		out_weights = redundancy_weights(
			idxs,
			use_zwt,
			use_dcurv,
			zwt_min,
			zwt_max,
			dcurv_min,
			dcurv_max);

		// === Class-balance 50/50 (SPRINGS vs OTHERS), exact proximity (EPS=1e-2) ===
		float Wspr = 0.0f, Woth = 0.0f;
		std::vector<uint8_t> is_spring(idxs.size(), 0);
		for (size_t k = 0; k < idxs.size(); ++k) {
			const int id = idxs[k];
			const bool sflag =
				is_spring_observation[static_cast<std::size_t>(id)] != 0u;

			is_spring[k] = sflag ? 1u : 0u;
			if (sflag) Wspr += out_weights[k];
			else Woth += out_weights[k];
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

	auto predict_with_beta = [&](const int node_id,
		const std::vector<float>& beta_local,
		const bool use_zwt_local,
		const bool use_dcurv_local,
		const float zwt_min_local,
		const float zwt_max_local,
		const float dcurv_min_local,
		const float dcurv_max_local) -> float
	{
		float prediction = beta_local[0];
		int coefficient_index = 1;

		if (use_zwt_local) {
			const float normalized_zwt =
				(zwt[node_id] - zwt_min_local) /
				std::max(1e-12f, zwt_max_local - zwt_min_local);
			prediction += beta_local[coefficient_index++] * normalized_zwt;
		}

		if (use_dcurv_local) {
			const float normalized_dcurv =
				(dcurv[node_id] - dcurv_min_local) /
				std::max(1e-12f, dcurv_max_local - dcurv_min_local);
			prediction += beta_local[coefficient_index++] * normalized_dcurv;
		}

		return prediction;
	};

	auto log_regression_state = [&](
		const char* stage,
		const std::vector<int>& indices,
		const std::vector<float>& beta_local,
		const std::vector<float>& weights_local,
		const bool active_zwt,
		const bool active_dcurv,
		const float zwt_min_local,
		const float zwt_max_local,
		const float dcurv_min_local,
		const float dcurv_max_local)
	{
		if (!K_EXT_DRIFT_ENABLE_DIAGNOSTIC_LOGS || indices.empty()) {
			return;
		}

		float y_min = std::numeric_limits<float>::infinity();
		float y_max = -std::numeric_limits<float>::infinity();
		double y_sum = 0.0;
		double weighted_y_sum = 0.0;
		double weight_sum = 0.0;

		int n_inlets = 0;
		int n_outlets = 0;
		int n_waypoints = 0;

		for (std::size_t k = 0; k < indices.size(); ++k) {
			const int id = indices[k];
			const float y = eq_radius_values[id];

			y_min = std::min(y_min, y);
			y_max = std::max(y_max, y);
			y_sum += static_cast<double>(y);

			if (k < weights_local.size()) {
				const double w =
					static_cast<double>(weights_local[k]);
				weighted_y_sum += w * static_cast<double>(y);
				weight_sum += w;
			}

			switch (conditioning_role_at(id)) {
			case ConditioningDataRole::Inlet:
				++n_inlets;
				break;
			case ConditioningDataRole::Outlet:
				++n_outlets;
				break;
			case ConditioningDataRole::Waypoint:
				++n_waypoints;
				break;
			default:
				break;
			}
		}

		std::cout
			<< "[geostats][external_drift][diag] " << stage
			<< " fit: n=" << indices.size()
			<< " (inlets=" << n_inlets
			<< ", outlets=" << n_outlets
			<< ", waypoints=" << n_waypoints << ")"
			<< "; Y[min/mean/max]=["
			<< y_min << ", "
			<< y_sum / static_cast<double>(indices.size()) << ", "
			<< y_max << "]";

		if (weight_sum > 0.0) {
			std::cout
				<< "; weighted Y mean="
				<< weighted_y_sum / weight_sum;
		}

		if (!beta_local.empty()) {
			std::cout << "; beta0=" << beta_local[0];

			int beta_index = 1;
			if (active_zwt) {
				std::cout
					<< ", beta_zwt="
					<< beta_local[beta_index++];
			}
			if (active_dcurv) {
				std::cout
					<< ", beta_dcurv="
					<< beta_local[beta_index++];
			}
		}

		if (active_zwt) {
			const auto full_bounds =
				std::minmax_element(zwt.begin(), zwt.end());

			std::cout
				<< "; zwt calib=[" << zwt_min_local
				<< ", " << zwt_max_local
				<< "], all=[" << *full_bounds.first
				<< ", " << *full_bounds.second << "]";
		}

		if (active_dcurv) {
			const auto full_bounds =
				std::minmax_element(dcurv.begin(), dcurv.end());

			std::cout
				<< "; dcurv calib=[" << dcurv_min_local
				<< ", " << dcurv_max_local
				<< "], all=[" << *full_bounds.first
				<< ", " << *full_bounds.second << "]";
		}

		std::cout << std::endl;
	};

	// ---- 2) Initial fit --------------------------------------------------------
	float zwt_min = 0.0f, zwt_max = 1.0f, dcurv_min = 0.0f, dcurv_max = 1.0f;
	std::vector<float> beta, weights_obs;
	if (!fit_on_subset(valid_indices, use_drift_zwt, use_drift_curv,
		beta, weights_obs, zwt_min, zwt_max, dcurv_min, dcurv_max)) {
		weights_out.assign(N, 0.0f);
		return drift;
	}

	log_regression_state(
		"initial",
		valid_indices,
		beta,
		weights_obs,
		use_drift_zwt,
		use_drift_curv,
		zwt_min,
		zwt_max,
		dcurv_min,
		dcurv_max
	);

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

	const bool drift_flags_changed =
		(drift_valid_zwt != use_drift_zwt || drift_valid_dcurv != use_drift_curv);

	if (K_EXT_DRIFT_ENABLE_DIAGNOSTIC_LOGS &&
		drift_flags_changed) {

		std::cout
			<< "[geostats][external_drift][diag] Geological sign check changed "
			<< "the active predictors: zwt="
			<< (drift_valid_zwt ? "kept" : "rejected")
			<< ", dcurv="
			<< (drift_valid_dcurv ? "kept" : "rejected")
			<< "."
			<< std::endl;
	}

	// Refit after a geological-sign rejection so that the retained coefficients
	// and normalization ranges are estimated with exactly the active predictors.
	if (drift_flags_changed && (drift_valid_zwt || drift_valid_dcurv)) {
		std::vector<float> beta_reduced;
		std::vector<float> weights_reduced;
		float zwt_min_reduced = zwt_min;
		float zwt_max_reduced = zwt_max;
		float dcurv_min_reduced = dcurv_min;
		float dcurv_max_reduced = dcurv_max;

		if (!fit_on_subset(
			valid_indices,
			drift_valid_zwt,
			drift_valid_dcurv,
			beta_reduced,
			weights_reduced,
			zwt_min_reduced,
			zwt_max_reduced,
			dcurv_min_reduced,
			dcurv_max_reduced))
		{
			weights_out.assign(N, 0.0f);
			return drift;
		}

		beta = beta_reduced;
		weights_obs = weights_reduced;
		zwt_min = zwt_min_reduced;
		zwt_max = zwt_max_reduced;
		dcurv_min = dcurv_min_reduced;
		dcurv_max = dcurv_max_reduced;
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
	if (K_EXT_DRIFT_ENABLE_MAD_TRIMMING) {
		int n_inlets = 0;
		int n_outlets = 0;
		for (const int id : valid_indices) {
			const ConditioningDataRole role = conditioning_role_at(id);
			if (role == ConditioningDataRole::Inlet) ++n_inlets;
			else if (role == ConditioningDataRole::Outlet) ++n_outlets;
		}

		const int n_var_cur =
			1 + (drift_valid_zwt ? 1 : 0) + (drift_valid_dcurv ? 1 : 0);
		std::vector<float> residuals_cur;
		residuals_cur.reserve(n_obs);
		std::vector<uint8_t> residual_is_testable;
		residual_is_testable.reserve(n_obs);

		// Cache local redundancy neighborhoods once. For ordinary omissions this
		// preserves the exact leave-one-out weighting while avoiding a complete
		// O(M^2) density reconstruction for each observation. Unique predictor
		// extrema deliberately fall back to the original complete subset fit.
		FastLeaveOneOutRegressionCache fast_loo_cache(
			valid_indices,
			zwt,
			dcurv,
			eq_radius_values,
			is_spring_observation,
			drift_valid_zwt,
			drift_valid_dcurv,
			zwt_min,
			zwt_max,
			dcurv_min,
			dcurv_max);

		for (int observation_position = 0; observation_position < n_obs; ++observation_position) {
			const int id = valid_indices[static_cast<std::size_t>(observation_position)];
			const ConditioningDataRole role = conditioning_role_at(id);
			bool loo_allowed = true;
			if (role == ConditioningDataRole::Inlet && n_inlets <= 1) {
				loo_allowed = false;
			}
			if (role == ConditioningDataRole::Outlet && n_outlets <= 1) {
				loo_allowed = false;
			}
			if (n_obs - 1 < n_var_cur) {
				loo_allowed = false;
			}

			float residual = 0.0f;
			bool used_loo = false;
			if (loo_allowed) {
				std::vector<float> beta_loo;
				float zwt_min_loo = zwt_min;
				float zwt_max_loo = zwt_max;
				float dcurv_min_loo = dcurv_min;
				float dcurv_max_loo = dcurv_max;
				bool fit_success = false;

				if (fast_loo_cache.can_fit_excluding(observation_position)) {
					fit_success = fast_loo_cache.fit_excluding(
						observation_position, beta_loo);
				}
				else {
					std::vector<int> loo_indices;
					loo_indices.reserve(valid_indices.size() - 1);
					for (const int other_id : valid_indices) {
						if (other_id != id) loo_indices.push_back(other_id);
					}
					std::vector<float> weights_loo;
					fit_success = fit_on_subset(
						loo_indices,
						drift_valid_zwt,
						drift_valid_dcurv,
						beta_loo,
						weights_loo,
						zwt_min_loo,
						zwt_max_loo,
						dcurv_min_loo,
						dcurv_max_loo);
				}

				if (fit_success) {
					const float prediction_loo = predict_with_beta(
						id,
						beta_loo,
						drift_valid_zwt,
						drift_valid_dcurv,
						zwt_min_loo,
						zwt_max_loo,
						dcurv_min_loo,
						dcurv_max_loo);
					residual = eq_radius_values[id] - prediction_loo;
					used_loo = true;
				}
			}

			if (!used_loo) {
				const float prediction_full = predict_with_beta(
					id,
					beta,
					drift_valid_zwt,
					drift_valid_dcurv,
					zwt_min,
					zwt_max,
					dcurv_min,
					dcurv_max);
				residual = eq_radius_values[id] - prediction_full;
			}

			residuals_cur.push_back(residual);
			residual_is_testable.push_back(used_loo ? 1u : 0u);
		}
		std::vector<float> absolute_residuals = residuals_cur;
		for (float& value : absolute_residuals) value = std::abs(value);
		std::nth_element(
			absolute_residuals.begin(),
			absolute_residuals.begin() + absolute_residuals.size() / 2,
			absolute_residuals.end());
		const float mad = absolute_residuals[absolute_residuals.size() / 2];
		const float sigma = std::max(1e-6f, 1.4826f * mad);

		// Two-sided central Gaussian compatibility interval of approximately 80%.
		//| Central fraction retained | Cutoff |
		//	| 80 % | 1.282 |
		//	| 85 % | 1.440 |
		//	| 90 % | 1.645 |
		//	| 95 % | 1.960 |
		//	| 97.5 % | 2.241 |
		//	| 98 % | 2.326 |
		//	| 99 % | 2.576 |
		//	| 99.5 % | 2.807 |
		//	| 99.9 % | 3.291 |
		const float C_CUTOFF = 2.241f;
		std::vector<int> survivors;
		survivors.reserve(n_obs);
		for (int observation = 0; observation < n_obs; ++observation) {
			const float standardized_residual =
				std::abs(residuals_cur[observation]) / sigma;
			if (!residual_is_testable[observation] ||
				standardized_residual <= C_CUTOFF) {
				survivors.push_back(valid_indices[observation]);
			}
		}

		if (K_EXT_DRIFT_ENABLE_DIAGNOSTIC_LOGS) {
			int trimmed_inlets = 0;
			int trimmed_outlets = 0;
			int trimmed_waypoints = 0;

			for (int observation = 0;
				observation < n_obs;
				++observation) {

				const float standardized_residual =
					std::abs(residuals_cur[observation]) / sigma;

				const bool trimmed =
					residual_is_testable[observation] &&
					standardized_residual > C_CUTOFF;

				if (!trimmed) {
					continue;
				}

				switch (conditioning_role_at(
					valid_indices[observation])) {

				case ConditioningDataRole::Inlet:
					++trimmed_inlets;
					break;
				case ConditioningDataRole::Outlet:
					++trimmed_outlets;
					break;
				case ConditioningDataRole::Waypoint:
					++trimmed_waypoints;
					break;
				default:
					break;
				}
			}

			std::cout
				<< "[geostats][external_drift][diag] LOO/MAD: MAD="
				<< mad
				<< ", robust sigma=" << sigma
				<< ", cutoff=" << C_CUTOFF
				<< "; retained=" << survivors.size()
				<< "/" << n_obs
				<< "; trimmed(inlets/outlets/waypoints)="
				<< trimmed_inlets << "/"
				<< trimmed_outlets << "/"
				<< trimmed_waypoints
				<< "."
				<< std::endl;
		}

		if (static_cast<int>(survivors.size()) >= n_var_cur &&
			static_cast<int>(survivors.size()) < n_obs) {
			std::vector<float> beta_refit;
			std::vector<float> weights_refit;
			float zwt_min_refit = zwt_min;
			float zwt_max_refit = zwt_max;
			float dcurv_min_refit = dcurv_min;
			float dcurv_max_refit = dcurv_max;

			if (fit_on_subset(
				survivors,
				drift_valid_zwt,
				drift_valid_dcurv,
				beta_refit,
				weights_refit,
				zwt_min_refit,
				zwt_max_refit,
				dcurv_min_refit,
				dcurv_max_refit))
			{
				beta = beta_refit;
				weights_obs = weights_refit;
				zwt_min = zwt_min_refit;
				zwt_max = zwt_max_refit;
				dcurv_min = dcurv_min_refit;
				dcurv_max = dcurv_max_refit;
				valid_indices = survivors;
			}
		}
	}

	log_regression_state(
		"final",
		valid_indices,
		beta,
		weights_obs,
		drift_valid_zwt,
		drift_valid_dcurv,
		zwt_min,
		zwt_max,
		dcurv_min,
		dcurv_max
	);

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

	if (K_EXT_DRIFT_ENABLE_DIAGNOSTIC_LOGS && !drift.empty()) {
		const auto bounds =
			std::minmax_element(drift.begin(), drift.end());

		const double mean =
			std::accumulate(
				drift.begin(),
				drift.end(),
				0.0
			) / static_cast<double>(drift.size());

		std::cout
			<< "[geostats][external_drift][diag] Final drift field: "
			<< "min=" << *bounds.first
			<< ", mean=" << mean
			<< ", max=" << *bounds.second
			<< "."
			<< std::endl;
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

	if (simulation_distribution == nullptr ||
		simulation_distribution->size() < 2u) {
		throw std::invalid_argument(
			"SGS requires a simulation distribution containing at least two values."
		);
	}

	std::unique_ptr<NormalScoreVariogramConverter> variogram_converter;
	if (!K_VARIOGRAM_PARAMETERS_ARE_ALREADY_GAUSSIAN) {
		variogram_converter = std::make_unique<NormalScoreVariogramConverter>(
			*simulation_distribution);
	}
	ScopedVariogramConverter scoped_variogram_converter(variogram_converter.get());

	// 1) Perform Normal Score Transform of initial distrib AND initial data vector
	std::vector<float> simulated_prop_gauss;
	if (!simulated_property.empty()) {
		simulated_prop_gauss = nst_data_with_nodata_fast_equivalent(simulated_property, *simulation_distribution, 0., 1.); // gaussianize the data values if any (do not change the no data values though)
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
			const Vector3& offending_point = curve->nodes.at(first_idx).p;

			throw std::runtime_error(
				"Sequential Gaussian Simulation aborted: conditioning value " +
				std::to_string(first_val) +
				" at skeleton node index " + std::to_string(first_idx) +
				" (zero-based; X=" + std::to_string(offending_point.x) +
				", Y=" + std::to_string(offending_point.y) +
				", Z=" + std::to_string(offending_point.z) +
				") lies outside the support [" +
				std::to_string(sim_min) + ", " +
				std::to_string(sim_max) +
				"] of the initial simulation distribution. " +
				"Total number of out-of-support conditioning data: " +
				std::to_string(n_outside) +
				". Please extend or change the initial simulation distribution accordingly."
			);
		}
	}

	//save_data(*simulation_distribution, "distrib.txt");
	//save_data(sim_distrib_gauss, "gauss_distrib.txt");

	// Cache branch membership once in global node-index order. The former code
	// rescanned the complete skeleton for every branch in both SGS passes. Keeping
	// each branch list ordered by node index preserves the exact input order passed
	// to `select_random_elements`, and therefore preserves the RNG sequence.
	std::vector<std::vector<int>> nodes_by_branch(curve->branch_sizes.size());
	for (int node_index = 0; node_index < static_cast<int>(curve->nodes.size()); ++node_index) {
		const int branch_id = curve->nodes[static_cast<std::size_t>(node_index)].branch_id;
		if (branch_id >= 0 && branch_id < static_cast<int>(nodes_by_branch.size())) {
			nodes_by_branch[static_cast<std::size_t>(branch_id)].push_back(node_index);
		}
	}

	// 2) Interbranch simulation with interbranch variogram

	// iterate on branches:
	for (int branch_id = 0; branch_id < curve->branch_sizes.size(); branch_id++) {

		// find number of nodes to simulate in given branch
		int nb_nodes_to_simulate = compute_prop_branch(curve->branch_sizes.at(branch_id), nb_points_interbranch, proportion_interbranch);

		std::vector<int> all_branch_nodes;
		all_branch_nodes.reserve(nodes_by_branch[static_cast<std::size_t>(branch_id)].size());
		for (int node_index : nodes_by_branch[static_cast<std::size_t>(branch_id)]) {
			if ((simulated_prop_gauss[static_cast<std::size_t>(node_index)] - (-99999)) < 1e-12) {
				all_branch_nodes.push_back(node_index);
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

			if (!neighbors_interbranch.empty() && std::isfinite(var_estimation)) { // conditional kriging succeeded
				var_estimation = (var_estimation < 0.) ? 0. : var_estimation; // due to numerical uncertainties, var estimation can sometimes be slightly negative.
				float simulated_val = generateNormalRandom(val_estimation, std::sqrt(var_estimation));
				simulated_prop_gauss.at(nodes_to_simulate.at(i)) = simulated_val;
			}
			else { // empty neighborhood or failed solve: use the single unconditional draw
				simulated_prop_gauss.at(nodes_to_simulate.at(i)) = val_estimation;
			}
		}
	}

	// 3) Intrabranch simulation with intrabranch variogram

	for (int branch_id = 0; branch_id < curve->branch_sizes.size(); branch_id++) {

		int nb_nodes_to_simulate = compute_prop_branch(curve->branch_sizes.at(branch_id), curve->branch_sizes.at(branch_id), 1.); // compute ALL remaining nodes

		std::vector<int> all_branch_nodes;
		all_branch_nodes.reserve(nodes_by_branch[static_cast<std::size_t>(branch_id)].size());
		for (int node_index : nodes_by_branch[static_cast<std::size_t>(branch_id)]) {
			if ((simulated_prop_gauss[static_cast<std::size_t>(node_index)] - (-99999)) < 1e-12) {
				all_branch_nodes.push_back(node_index);
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


			if (!neighbors_intrabranch.empty() && std::isfinite(var_estimation)) {
				var_estimation = (var_estimation < 0.) ? 0. : var_estimation;// due to numerical uncertainties, var estimation can sometimes be slightly negative.
				float simulated_val = generateNormalRandom(val_estimation, std::sqrt(var_estimation));
				simulated_prop_gauss.at(nodes_to_simulate.at(i)) = simulated_val;
			}
			else { // empty neighborhood or failed solve: use the single unconditional draw
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

		if (!neighbors_global.empty() && std::isfinite(var_estimation_global)) {
			var_estimation_global = (var_estimation_global < 0.) ? 0. : var_estimation_global;// due to numerical uncertainties, var estimation can sometimes be slightly negative.
			float simulated_val_global = generateNormalRandom(val_estimation_global, std::sqrt(var_estimation_global));
			simulated_prop_gauss.at(nodes_to_simulate_global.at(i)) = simulated_val_global;
		}
		else { // empty neighborhood or failed solve: use the single unconditional draw
			simulated_prop_gauss.at(nodes_to_simulate_global.at(i)) = val_estimation_global;
		}
	}

	// 5) Back-transform of generated distribution.

	// if we had data at the beginning, use it as basis for back transformation
	//if (!simulated_property.empty()) {
	//	simulated_property = back_transform(simulated_prop_gauss, simulated_property_copy);
	//}
	//else { // else, use the initial distribution instead
	simulated_property = back_transform_fast_equivalent(simulated_prop_gauss, *simulation_distribution);
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
	const std::vector<ConditioningDataRole>& conditioning_roles,
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
		conditioning_roles,
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
	// === 3) Compute conditioning residuals and estimate background residual scale ===
	//
	// Every hard conditioning datum remains in the residual vector so that SGS can
	// honor it locally, including observations excluded from the drift regression.
	// The background residual variance is estimated only from observations retained
	// in the final drift fit. `weights_output` is zero for observations excluded by
	// the active regression filters or robust trimming and positive for final-fit data.
	std::vector<float> residuals(simulated_property.size(), -99999.0f);
	double weighted_residual_squared_sum = 0.0;
	double residual_weight_sum = 0.0;
	float conditioning_residual_min = std::numeric_limits<float>::infinity();
	float conditioning_residual_max = -std::numeric_limits<float>::infinity();
	std::size_t conditioning_residual_count = 0;

	for (size_t i = 0; i < simulated_property.size(); ++i) {
		if (std::abs(simulated_property[i] - (-99999.0f)) <= 1e-12f) {
			continue;
		}

		const float residual = simulated_property[i] - drift[i];
		if (!std::isfinite(residual)) {
			throw std::runtime_error(
				"[geostats][external_drift] A non-finite conditioning residual was produced."
			);
		}

		residuals[i] = residual;
		conditioning_residual_min = std::min(conditioning_residual_min, residual);
		conditioning_residual_max = std::max(conditioning_residual_max, residual);
		++conditioning_residual_count;

		if (i < weights_output.size() && weights_output[i] > 0.0f) {
			const double weight = static_cast<double>(weights_output[i]);
			weighted_residual_squared_sum +=
				weight * static_cast<double>(residual) * static_cast<double>(residual);
			residual_weight_sum += weight;
		}
	}

	if (conditioning_residual_count == 0) {
		throw std::runtime_error(
			"[geostats][external_drift] Cannot construct the residual simulation: "
			"no valid conditioning observation is available."
		);
	}

	if (!std::isfinite(residual_weight_sum) || residual_weight_sum <= 0.0) {
		throw std::runtime_error(
			"[geostats][external_drift] Cannot estimate background residual variability: "
			"the final drift fit contains no positively weighted observation."
		);
	}

	// The final WLS fit includes an intercept, so its residual process is modeled
	// around zero. Use the same final regression weights to estimate the unexplained
	// background variability without allowing deliberately excluded hard anomalies
	// to inflate the residual variance throughout the network.
	const double sigma_residual_background = std::sqrt(
		weighted_residual_squared_sum / residual_weight_sum
	);

	// === 4) Build the background residual distribution and add rare tail anchors ===
	//
	// Preserve the standardized shape of the supplied total-property marginal while
	// replacing its spread by the residual variability unexplained by the final drift:
	//
	//     e_j = (Y_j - mean(Y)) * sigma_e / sigma_Y
	//
	// Hard conditioning residuals excluded from the drift fit must still be honored.
	// If they lie beyond the background support, append only the most extreme lower
	// and/or upper residual as tail anchors. This extends the numerical support without
	// letting the number of anomalous conditioning points define their probability mass.
	if (simulation_distribution == nullptr || simulation_distribution->empty()) {
		throw std::invalid_argument(
			"[geostats][external_drift] A non-empty simulation distribution is required "
			"to construct the residual distribution."
		);
	}

	const double distribution_mean = std::accumulate(
		simulation_distribution->begin(),
		simulation_distribution->end(),
		0.0
	) / static_cast<double>(simulation_distribution->size());

	double distribution_squared_deviation_sum = 0.0;
	for (const float value : *simulation_distribution) {
		if (!std::isfinite(value)) {
			throw std::runtime_error(
				"[geostats][external_drift] Cannot construct the residual distribution: "
				"the supplied simulation distribution contains a non-finite value."
			);
		}
		const double centered = static_cast<double>(value) - distribution_mean;
		distribution_squared_deviation_sum += centered * centered;
	}

	const double sigma_distribution = std::sqrt(
		distribution_squared_deviation_sum /
		static_cast<double>(simulation_distribution->size())
	);

	if (!std::isfinite(sigma_distribution) ||
		sigma_distribution <= std::numeric_limits<double>::epsilon()) {
		throw std::runtime_error(
			"[geostats][external_drift] Cannot construct the residual distribution: "
			"the supplied simulation distribution has zero or non-finite variance."
		);
	}

	if (!std::isfinite(sigma_residual_background)) {
		throw std::runtime_error(
			"[geostats][external_drift] Cannot construct the residual distribution: "
			"the background residual variance is non-finite."
		);
	}

	const double residual_scale =
		sigma_residual_background / sigma_distribution;

	std::vector<float> residual_sim_distribution;
	residual_sim_distribution.reserve(simulation_distribution->size() + 2u);
	for (const float value : *simulation_distribution) {
		const double centered = static_cast<double>(value) - distribution_mean;
		residual_sim_distribution.push_back(
			static_cast<float>(centered * residual_scale)
		);
	}

	auto background_bounds = std::minmax_element(
		residual_sim_distribution.begin(), residual_sim_distribution.end());
	const float background_min = *background_bounds.first;
	const float background_max = *background_bounds.second;

	if (K_EXT_DRIFT_ENABLE_DIAGNOSTIC_LOGS) {
		std::cout
			<< "[geostats][external_drift][diag] Residual SGS: "
			<< "sigma(property)=" << sigma_distribution
			<< ", sigma(background residual)="
			<< sigma_residual_background
			<< ", residual scale=" << residual_scale
			<< "; hard residual range=["
			<< conditioning_residual_min << ", "
			<< conditioning_residual_max << "]"
			<< "; background residual support=["
			<< background_min << ", "
			<< background_max << "]"
			<< "; tail anchors(lower/upper)="
			<< (conditioning_residual_min < background_min ? "yes" : "no")
			<< "/"
			<< (conditioning_residual_max > background_max ? "yes" : "no")
			<< "."
			<< std::endl;
	}

	if (conditioning_residual_min < background_min) {
		residual_sim_distribution.push_back(conditioning_residual_min);
	}
	if (conditioning_residual_max > background_max) {
		residual_sim_distribution.push_back(conditioning_residual_max);
	}

	std::sort(residual_sim_distribution.begin(), residual_sim_distribution.end());

	SGS3(
		curve,
		residuals,
		&residual_sim_distribution, // background residual distribution with rare hard-data tail anchors
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
