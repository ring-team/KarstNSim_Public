/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

***************************************************************/

#include "KarstNSim/randomgenerator.h"

std::mt19937 globalRng;
std::vector<std::uint32_t> globalSeed;

void initializeRng(const std::vector<std::uint32_t>& seed) {
	globalSeed = seed;
	globalRng = std::mt19937(globalSeed[0]); // Use the first element of the seed vector
}

float generateNormalRandom(float mean, float stddev) {
	std::normal_distribution<float> normalDistribution(mean, stddev);
	return normalDistribution(globalRng);
}

int generateRandomInt(int min, int max) {
	std::uniform_int_distribution<int> distribution(min, max);
	return distribution(globalRng);
}

float generateRandomFloat(float min, float max) {
	std::uniform_real_distribution<float> distribution(min, max);
	return distribution(globalRng);
}

bool generateRandomBool() {
	std::bernoulli_distribution distribution(0.5);
	return distribution(globalRng);
}

char generateRandomChar() {
	// For simplicity, generating a random ASCII character between 'a' and 'z'
	std::uniform_int_distribution<int> distribution('a', 'z');
	return static_cast<char>(distribution(globalRng));
}

std::size_t generateRandomIndex(std::size_t containerSize) {
	std::uniform_int_distribution<std::size_t> distribution(0, containerSize - 1);
	return distribution(globalRng);
}

std::vector<float> generateGaussianVector(std::size_t size, float mean, float stddev) {
	std::vector<float> result;
	result.reserve(size);

	std::normal_distribution<float> distribution(mean, stddev);
	for (std::size_t i = 0; i < size; ++i) {
		result.push_back(distribution(globalRng));
	}

	return result;
}

std::vector<int> generateRandomIntVector(std::size_t size, int minValue, int maxValue) {
	std::vector<int> result;
	result.reserve(size);

	std::uniform_int_distribution<int> distribution(minValue, maxValue);
	for (std::size_t i = 0; i < size; ++i) {
		result.push_back(distribution(globalRng));
	}

	return result;
}

std::vector<float> generateRandomFloatVector(std::size_t size, float minValue, float maxValue) {
	std::vector<float> result;
	result.reserve(size);

	std::uniform_real_distribution<float> distribution(minValue, maxValue);
	for (std::size_t i = 0; i < size; ++i) {
		result.push_back(distribution(globalRng));
	}

	return result;
}

// Function to compute the CDF of the standard normal distribution
float normal_cdf(float x) {
	return 0.5f * std::erfc(-x / std::sqrt(2.0f));
}

// Approximate inverse of the complementary error function (erfcinv) using Newton's method
float approximate_erfcinv(float x) {
	float guess = 0.5f;
	float tolerance = 1e-6f;
	int maxIterations = 100;

	for (int i = 0; i < maxIterations; ++i) {
		float error = std::erfc(guess) - x;
		if (std::abs(error) < tolerance)
			break;
		float derivative = static_cast<float>(-2.0f / std::sqrt(M_PI) * std::exp(-guess * guess));
		guess -= error / derivative;
	}

	return guess;
}

// Function to compute the inverse CDF of the standard normal distribution
float inverse_normal_cdf(float p) {
	return std::sqrt(2.0f) * approximate_erfcinv(2 * p);
}

// Function to perform Gaussian transform on a discrete distribution
std::vector<float> nst(const std::vector<float>& discreteDistribution, float mean, float stddev) {
	std::size_t size = discreteDistribution.size();
	std::vector<float> transformedValues(size);

	// Compute the ranks of the data
	std::vector<size_t> ranks(size);
	std::iota(ranks.begin(), ranks.end(), 0);
	std::sort(ranks.begin(), ranks.end(), [&](size_t i, size_t j) { return discreteDistribution[i] > discreteDistribution[j]; }); // Reverse order

	// Compute the quantiles of the original distribution
	std::vector<float> quantiles(size);
	for (std::size_t i = 0; i < size; ++i) {
		quantiles[ranks[i]] = (i + 0.5f) / size;
	}

	// Adjust quantiles for elements with the same value (duplicates) in the original discrete distribution
	for (std::size_t i = 1; i < size; ++i) {
		if (discreteDistribution[ranks[i]] == discreteDistribution[ranks[i - 1]]) {
			quantiles[ranks[i]] = quantiles[ranks[i - 1]]; // Assign same quantile
		}
	}

	// Compute Gaussian random values corresponding to the quantiles
	for (std::size_t i = 0; i < size; ++i) {
		transformedValues[i] = mean + stddev * inverse_normal_cdf(quantiles[i]);
	}

	return transformedValues;
}

// Function to perform Gaussian transform on a data vector, based on a discrete distribution (data values HAVE to be in the discrete distribution (otherwise
// the found quantile will be >1 or <0, and the corresponding gaussian value +inf or -inf respectively)
std::vector<float> nst_data_with_nodata(const std::vector<float>& data_vector, const std::vector<float>& discreteDistribution, float mean, float stddev) {
	std::size_t size = data_vector.size();
	std::vector<float> transformedValues(size);

	// Create a copy of the input distribution excluding "no data" values
	std::vector<float> filteredDistribution;
	for (const auto& value : discreteDistribution) {
		if (std::abs(value - (-99999)) > 1e-12) {
			filteredDistribution.push_back(value);
		}
	}

	// Compute the ranks of the filtered data
	std::vector<size_t> ranks(filteredDistribution.size());
	std::iota(ranks.begin(), ranks.end(), 0);
	std::sort(ranks.begin(), ranks.end(), [&](size_t i, size_t j) { return filteredDistribution[i] > filteredDistribution[j]; });

	// Compute the quantiles of the original distribution
	std::vector<float> quantiles(filteredDistribution.size());
	for (std::size_t i = 0; i < filteredDistribution.size(); ++i) {
		quantiles[ranks[i]] = (i + 0.5f) / filteredDistribution.size();
	}

	// Adjust quantiles for elements with the same value (duplicates) in the original discrete distribution
	for (std::size_t i = 1; i < filteredDistribution.size(); ++i) {
		if (filteredDistribution[ranks[i]] == filteredDistribution[ranks[i - 1]]) {
			quantiles[ranks[i]] = quantiles[ranks[i - 1]]; // Assign same quantile
		}
	}

	// Map the quantiles to the corresponding Gaussian random values
	for (std::size_t i = 0; i < size; ++i) {
		if (std::abs(data_vector[i] - (-99999)) > 1e-12) {
			// find the corresponding quantile in the distribution to the current data
			float data = data_vector[i];
			float q;
			// find the two closest quantiles in the initial discrete distribution, and interpolate between them.
			// make sure that the closest are one greater, one smaller, otherwise we extrapolate with linear model instead based on just the closest quantile value
			bool hasTwoNeighbors;
			//float tested_quantile = quantiles[j];
			std::pair<std::pair<float, size_t>, std::pair<float, size_t>> closest_vals = find_closest_values(filteredDistribution, data, hasTwoNeighbors);

			if (hasTwoNeighbors) { // interpolation case
				float q1 = quantiles[ranks[closest_vals.first.second]];
				float q2 = quantiles[ranks[closest_vals.second.second]];
				q = interpolate(q1, closest_vals.first.first, q2, closest_vals.second.first, data);
			}
			else { // extrapolate (or already existing value) case
				float q1 = quantiles[ranks[closest_vals.first.second]];
				float q2 = quantiles[ranks[closest_vals.second.second]];
				if (closest_vals.first.first < data) { // case where tested quantile is larger than any value in the initial distrib
					q = interpolate(q1, closest_vals.first.first, q2, closest_vals.second.first, data);
				}
				else if (closest_vals.first.first > data) { // case where tested quantile is smaller than any value in the initial distrib
					q = interpolate(q1, closest_vals.first.first, q2, closest_vals.second.first, data);
				}
				else { // case where closest_vals.first.first == data (value already in inital distribution)
					q = q1;
				}
			}
			transformedValues[i] = mean + stddev * inverse_normal_cdf(1 - q);
		}
		else {
			transformedValues[i] = -99999;
		}
	}

	return transformedValues;
}

// Function to perform Normal Score Transform on a discrete distribution with no data values in it (that shouldnt be transformed)
std::vector<float> nst_with_nodata(const std::vector<float>& discreteDistribution, float mean, float stddev) {
	std::size_t size = discreteDistribution.size();
	std::vector<float> transformedValues(size);

	// Create a copy of the input distribution excluding "no data" values
	std::vector<float> filteredDistribution;
	for (const auto& value : discreteDistribution) {
		if (std::abs(value - (-99999)) > 1e-12) {
			filteredDistribution.push_back(value);
		}
	}

	// Compute the ranks of the filtered data
	std::vector<size_t> ranks(filteredDistribution.size());
	std::iota(ranks.begin(), ranks.end(), 0);
	std::sort(ranks.begin(), ranks.end(), [&](size_t i, size_t j) { return filteredDistribution[i] > filteredDistribution[j]; });

	// Compute the quantiles of the original distribution
	std::vector<float> quantiles(filteredDistribution.size());
	for (std::size_t i = 0; i < filteredDistribution.size(); ++i) {
		quantiles[ranks[i]] = (i + 0.5f) / filteredDistribution.size();
	}

	// Adjust quantiles for elements with the same value (duplicates) in the original discrete distribution
	for (std::size_t i = 1; i < filteredDistribution.size(); ++i) {
		if (filteredDistribution[ranks[i]] == filteredDistribution[ranks[i - 1]]) {
			quantiles[ranks[i]] = quantiles[ranks[i - 1]]; // Assign same quantile
		}
	}

	// Map the quantiles to the corresponding Gaussian random values
	for (std::size_t i = 0, j = 0; i < size; ++i) {
		if (std::abs(discreteDistribution[i] - (-99999)) > 1e-12) {
			transformedValues[i] = mean + stddev * inverse_normal_cdf(quantiles[j++]);
		}
		else {
			transformedValues[i] = -99999;
		}
	}

	return transformedValues;
}

std::pair<std::pair<float, size_t>, std::pair<float, size_t>> find_closest_values(const std::vector<float>& vec, float value, bool& hasTwoNeighbors) {
	std::vector<std::pair<float, size_t>> indexedVec;
	for (size_t i = 0; i < vec.size(); ++i) {
		indexedVec.emplace_back(vec[i], i); // Store each value along with its index
	}

	std::sort(indexedVec.begin(), indexedVec.end(), [&](const auto& lhs, const auto& rhs) {
		float lhs_diff = std::abs(lhs.first - value);
		float rhs_diff = std::abs(rhs.first - value);
		if (lhs_diff == rhs_diff) {
			return lhs.first < rhs.first; // If equidistant, prioritize the smaller value
		}
		return lhs_diff < rhs_diff; // Sort based on absolute difference to value
	});

	if (indexedVec.empty()) {
		hasTwoNeighbors = false;
		return { {0.0, 0}, {0.0, 0} }; // Return default values if vector is empty
	}

	// Check if the closest neighbor is equal to the input value
	if (std::abs(indexedVec[0].first - value) < 1e-4) {
		hasTwoNeighbors = false;
	}
	else {
		// Check if both neighbors are greater than or smaller than the input value
		if ((indexedVec[0].first > value && indexedVec[1].first > value) || (indexedVec[0].first < value && indexedVec[1].first < value)) {
			hasTwoNeighbors = false;
		}
		else {
			hasTwoNeighbors = true;
		}
	}

	// If the first element is closest to value, return it and the second closest
	return { {indexedVec[0].first, indexedVec[0].second}, {indexedVec[1].first, indexedVec[1].second} };
}


// Function to perform linear interpolation
float interpolate(float x1, float y1, float x2, float y2, float y) {
	// Check if y1 and y2 are equal
	if (y1 == y2)
		return (x1 + x2) / 2.0f; // Return the average of x1 and x2

	// Calculate x using linear interpolation formula
	return x1 + ((x2 - x1) * (y - y1)) / (y2 - y1);
}

// Function to perform the back transform from Gaussian distribution to discrete distribution (important precision: I chose to interpolate (and extrapolate) values linearly!!)
std::vector<float> back_transform(const std::vector<float>& gaussianDistribution, const std::vector<float>& discreteDistribution) {
	std::size_t size = gaussianDistribution.size();
	std::vector<float> backTransformedValues(size);

	// 1) Find quantiles on gaussian distribution

	// Create a copy of the input gaussian distribution excluding "no data" values
	std::vector<float> filteredDistribution;
	for (const auto& value : gaussianDistribution) {
		if (std::abs(value - (-99999)) > 1e-12) {
			filteredDistribution.push_back(value);
		}
	}
	// Compute the quantiles of the input gaussian distribution
	std::vector<float> quantiles(filteredDistribution.size());
	for (std::size_t i = 0; i < filteredDistribution.size(); ++i) {
		quantiles[i] = normal_cdf(filteredDistribution[i]);
	}

	// 2) Find quantiles on initial discrete distribution

	// Create a copy of the input discrete distribution excluding "no data" values
	std::vector<float> filteredInitialDistribution;
	for (const auto& value : discreteDistribution) {
		if (std::abs(value - (-99999)) > 1e-12) {
			filteredInitialDistribution.push_back(value);
		}
	}
	// Compute the ranks of the filtered data
	std::vector<size_t> ranksInitial(filteredInitialDistribution.size());
	std::iota(ranksInitial.begin(), ranksInitial.end(), 0);
	std::sort(ranksInitial.begin(), ranksInitial.end(), [&](size_t i, size_t j) { return filteredInitialDistribution[i] < filteredInitialDistribution[j]; });
	// Compute the quantiles of the input discrete distribution
	std::vector<float> quantilesInitial(filteredInitialDistribution.size());
	for (std::size_t i = 0; i < filteredInitialDistribution.size(); ++i) {
		quantilesInitial[ranksInitial[i]] = (i + 0.5f) / filteredInitialDistribution.size();
	}

	// 3) Map the quantiles to the corresponding initial distribution, with interpolation or extrapolation
	for (std::size_t i = 0, j = 0; i < size; ++i) {
		if (std::abs(gaussianDistribution[i] - (-99999)) > 1e-12) {
			// find the two closest quantiles in the initial discrete distribution, and interpolate between them.
			// make sure that the closest are one greater, one smaller, otherwise we extrapolate with linear model instead based on just the closest quantile value
			bool hasTwoNeighbors;
			float tested_quantile = quantiles[j];
			std::pair<std::pair<float, size_t>, std::pair<float, size_t>> closest_quantiles = find_closest_values(quantilesInitial, quantiles[j++], hasTwoNeighbors);

			if (hasTwoNeighbors) { // interpolation case
				float val1 = filteredInitialDistribution.at(closest_quantiles.first.second);
				float val2 = filteredInitialDistribution.at(closest_quantiles.second.second);
				backTransformedValues[i] = interpolate(val1, closest_quantiles.first.first, val2, closest_quantiles.second.first, tested_quantile);
			}
			else { // extrapolate (or already existing value) case
				float val1 = filteredInitialDistribution.at(closest_quantiles.first.second);
				float val2 = filteredInitialDistribution.at(closest_quantiles.second.second);
				if (closest_quantiles.first.first < tested_quantile) { // case where tested quantile is larger than any value in the initial distrib
					backTransformedValues[i] = interpolate(val1, closest_quantiles.first.first, val2, closest_quantiles.second.first, tested_quantile);
				}
				else if (closest_quantiles.first.first > tested_quantile) { // case where tested quantile is smaller than any value in the initial distrib
					backTransformedValues[i] = interpolate(val1, closest_quantiles.first.first, val2, closest_quantiles.second.first, tested_quantile);
				}
				else { // case where closest_quantiles.first.first == tested_quantile (value already in inital distribution)
					backTransformedValues[i] = val1;
				}
			}
			// if get there, we defined the new value
		}
		else {
			backTransformedValues[i] = -99999;
		}
	}

	return backTransformedValues;
}

// Function to select n elements out of k in a vector at random
float select_random_element(const std::vector<float>& vec) {
	float result;

	// Shuffle the input vector
	std::vector<float> shuffledVec = vec;
	std::shuffle(shuffledVec.begin(), shuffledVec.end(), globalRng);

	// Select the first n elements from the shuffled vector
	result = shuffledVec[0];

	return result;
}

// Function to select n elements out of k in a vector at random
std::vector<int> select_random_elements(const std::vector<int>& vec, int n) {

	std::vector<int> result;

	// Check if n is greater than the size of the vector
	if (n > int(vec.size())) {
		n = int(vec.size());
	}

	// Shuffle the input vector
	std::vector<int> shuffledVec = vec;
	std::shuffle(shuffledVec.begin(), shuffledVec.end(), globalRng);

	// Select the first n elements from the shuffled vector
	for (int i = 0; i < n; ++i) {
		result.push_back(shuffledVec[i]);
	}

	return result;
}

// Function to compute the CDF of a discrete distribution
std::vector<float> computeCDF(const std::vector<float>& probabilities) {
	std::vector<float> cdf(probabilities.size());
	std::partial_sum(probabilities.begin(), probabilities.end(), cdf.begin());
	return cdf;
}

// Function to generate a random value from a discrete distribution using the CDF
int generateRandomDiscrete(const std::vector<float>& probabilities) {

	// Generate a random number between 0 and 1
	float rand_num = generateRandomFloat(0., 1.);

	// Compute the CDF of the discrete distribution
	std::vector<float> cdf = computeCDF(probabilities);

	// Find the first index in the CDF where the value exceeds the random number
	for (int i = 0; i < cdf.size(); ++i) {
		if (rand_num < cdf[i]) {
			return i;
		}
	}
	return int(cdf.size()) - 1; // In case rand_num is very close to 1
}

// Function to normalize a distribution into a probability distribution
std::vector<float> normalizeDistribution(const std::vector<float>& values) {
	// Find the minimum and maximum values in the distribution
	float min_value = *std::min_element(values.begin(), values.end());
	float max_value = *std::max_element(values.begin(), values.end());

	// Normalize the values to range between 0 and 1
	std::vector<float> normalized;
	for (float value : values) {
		float normalized_value = (value - min_value) / (max_value - min_value);
		normalized.push_back(normalized_value);
	}

	return normalized;
}