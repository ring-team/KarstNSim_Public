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
	for (int i = 0; i < N; ++i) {
		if (i < (int)eq_radius_values.size()) {
			const float r = eq_radius_values[i];
			const bool has_data = (std::abs(r - (-99999.0f)) > 1e-12f);
			const bool pass_radius_cap =
				(!K_EXT_DRIFT_ENABLE_RADIUS_CAP || r <= K_RADIUS_MAX_FOR_REGRESSION);
			if (has_data && pass_radius_cap) {
				valid_indices.push_back(i);
			}
		}
	}

	const int n_obs = static_cast<int>(valid_indices.size());
	if (n_obs == 0 || (!use_drift_zwt && !use_drift_curv)) {
		weights_out.assign(N, 0.0f);
		return drift;
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
			for (int k = 0; k < M; ++k) {
				std::vector<std::pair<float, int>> local_dist;
				local_dist.reserve(M);

				for (int s = 0; s < M; ++s) {
					const float dz = z01[s] - z01[k];
					const float dd = d01[s] - d01[k];
					local_dist.emplace_back(std::sqrt(dz * dz + dd * dd), s);
				}

				std::sort(local_dist.begin(), local_dist.end());

				float acc = 0.0f, norm = 0.0f;
				for (int s = 0; s < local_count; ++s) {
					const float distance = local_dist[s].first;
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
		const float SPR_EPS = 1e-2f;
		float Wspr = 0.0f, Woth = 0.0f;
		std::vector<uint8_t> is_spring(idxs.size(), 0);
		for (size_t k = 0; k < idxs.size(); ++k) {
			const int id = idxs[k];
			const Vector3& p = curve->nodes[id].p;

			bool sflag = false;
			for (size_t si = 0; si < springs_xyz.size(); ++si) {
				const float dist = KarstNSim::magnitude(p - springs_xyz[si]);
				if (dist <= SPR_EPS) { sflag = true; break; }
			}

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

	const bool drift_flags_changed =
		(drift_valid_zwt != use_drift_zwt || drift_valid_dcurv != use_drift_curv);

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

		for (const int id : valid_indices) {
			const ConditioningDataRole role = conditioning_role_at(id);
			bool loo_allowed = true;
			if (role == ConditioningDataRole::Inlet && n_inlets <= 1) {
				loo_allowed = false;
			}
			if (role == ConditioningDataRole::Outlet && n_outlets <= 1) {
				loo_allowed = false;
			}

			std::vector<int> loo_indices;
			loo_indices.reserve(
				valid_indices.empty() ? 0 : valid_indices.size() - 1);
			for (const int other_id : valid_indices) {
				if (other_id != id) loo_indices.push_back(other_id);
			}
			if (static_cast<int>(loo_indices.size()) < n_var_cur) {
				loo_allowed = false;
			}

			float residual = 0.0f;
			bool used_loo = false;
			if (loo_allowed) {
				std::vector<float> beta_loo;
				std::vector<float> weights_loo;
				float zwt_min_loo = zwt_min;
				float zwt_max_loo = zwt_max;
				float dcurv_min_loo = dcurv_min;
				float dcurv_max_loo = dcurv_max;

				if (fit_on_subset(
					loo_indices,
					drift_valid_zwt,
					drift_valid_dcurv,
					beta_loo,
					weights_loo,
					zwt_min_loo,
					zwt_max_loo,
					dcurv_min_loo,
					dcurv_max_loo))
				{
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
		const float C_CUTOFF = 1.28f;
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
