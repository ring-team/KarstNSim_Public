/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

***************************************************************/

#pragma once

/**
@file randomgenerator.h
@brief Functions allowing generation of seeds, random numbers and random processes.
@author Augustin GOUY
**/

#include <random>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Declare the random number generator and seed globally
extern std::mt19937 globalRng;
extern std::vector<std::uint32_t> globalSeed;

/**
 * @brief Initializes the random number generator with a given seed.
 * @param seed Vector of integers used as the seed.
 */
void initializeRng(const std::vector<std::uint32_t>& seed);

/**
 * @brief Generates a random number from a normal (Gaussian) distribution.
 * @param mean Mean of the distribution.
 * @param stddev Standard deviation of the distribution.
 * @return Random number from the specified normal distribution.
 */
float generateNormalRandom(float mean, float stddev);

/**
 * @brief Generates a random integer within a given range.
 * @param min Minimum value of the range (inclusive).
 * @param max Maximum value of the range (inclusive).
 * @return Random integer within the specified range.
 */
int generateRandomInt(int min, int max);

/**
 * @brief Generates a random float within a given range.
 * @param min Minimum value of the range (inclusive).
 * @param max Maximum value of the range (inclusive).
 * @return Random float within the specified range.
 */
float generateRandomFloat(float min, float max);

/**
 * @brief Generates a random boolean value (true or false).
 * @return Random boolean value.
 */
bool generateRandomBool();

/**
 * @brief Generates a random character.
 * @return Random character from the standard ASCII set.
 */
char generateRandomChar();

/**
 * @brief Generates a random index within a given container size.
 * @param containerSize Size of the container.
 * @return Random index in the range [0, containerSize-1].
 */
std::size_t generateRandomIndex(std::size_t containerSize);

/**
 * @brief Shuffles a container using the global random number generator.
 * @tparam ForwardIt Iterator type for the container.
 * @param first Iterator to the beginning of the container.
 * @param last Iterator to the end of the container.
 */
template <typename ForwardIt>
void shuffleContainer(ForwardIt first, ForwardIt last) {
	std::shuffle(first, last, globalRng);
}

/**
 * @brief Generates a vector of random numbers following a Gaussian distribution.
 * @param size Size of the vector.
 * @param mean Mean of the Gaussian distribution.
 * @param stddev Standard deviation of the Gaussian distribution.
 * @return Vector of random numbers.
 */
std::vector<float> generateGaussianVector(std::size_t size, float mean, float stddev);

/**
 * @brief Generates a vector of random integers within a given range.
 * @param size Size of the vector.
 * @param minValue Minimum value of the range.
 * @param maxValue Maximum value of the range.
 * @return Vector of random integers.
 */
std::vector<int> generateRandomIntVector(std::size_t size, int minValue, int maxValue);

/**
 * @brief Generates a vector of random floats within a given range.
 * @param size Size of the vector.
 * @param minValue Minimum value of the range.
 * @param maxValue Maximum value of the range.
 * @return Vector of random floats.
 */
std::vector<float> generateRandomFloatVector(std::size_t size, float minValue, float maxValue);

/**
 * @brief Transforms a distribution into a Gaussian distribution N(mean, stddev^2).
 * @param discreteDistribution Input discrete distribution.
 * @param mean Mean of the target Gaussian distribution.
 * @param stddev Standard deviation of the target Gaussian distribution.
 * @return Transformed Gaussian distribution.
 */
std::vector<float> nst(const std::vector<float>& discreteDistribution, float mean, float stddev);

/**
 * @brief Transforms a distribution to a Gaussian distribution, keeping no-data values (-99999) unchanged.
 * @param discreteDistribution Input discrete distribution.
 * @param mean Mean of the target Gaussian distribution.
 * @param stddev Standard deviation of the target Gaussian distribution.
 * @return Transformed Gaussian distribution.
 */
std::vector<float> nst_with_nodata(const std::vector<float>& discreteDistribution, float mean, float stddev);

/**
 * @brief Function to perform Gaussian transform on a data vector, based on a discrete distribution (data values HAVE to be in the discrete distribution, otherwise the found quantile will be >1 or <0, and the corresponding gaussian value +inf or -inf respectively)
 * @param data_vector Data vector to transform.
 * @param discreteDistribution Discrete distribution used for transformation.
 * @param mean Mean of the Gaussian distribution.
 * @param stddev Standard deviation of the Gaussian distribution.
 * @return Transformed data vector.
 */
std::vector<float> nst_data_with_nodata(const std::vector<float>& data_vector, const std::vector<float>& discreteDistribution, float mean, float stddev);

/**
 * @brief Finds the two closest values and their indices in a vector to a given value.
 * @param vec Input vector.
 * @param value Value to find neighbors for.
 * @param hasTwoNeighbors Boolean indicating if two neighbors were found.
 * @return Pair of closest values and their indices.
 */
std::pair<std::pair<float, size_t>, std::pair<float, size_t>> find_closest_values(const std::vector<float>& vec, float value, bool& hasTwoNeighbors);

/**
 * @brief Interpolates a value given two surrounding points.
 * @param x Input x value.
 * @param x0 x-coordinate of the first point.
 * @param x1 x-coordinate of the second point.
 * @param y0 y-coordinate of the first point.
 * @param y1 y-coordinate of the second point.
 * @return Interpolated y value.
 */
float interpolate(float x, float x0, float x1, float y0, float y1);

/**
 * @brief Back-transforms a Gaussian distribution to its original distribution.
 * @param gaussianDistribution Gaussian distribution to back-transform.
 * @param discreteDistribution Original discrete distribution.
 * @return Back-transformed distribution.
 */
std::vector<float> back_transform(const std::vector<float>& gaussianDistribution, const std::vector<float>& discreteDistribution);

/**
 * @brief Selects a random element from a vector.
 * @param vec Input vector.
 * @return Randomly selected element.
 */
float select_random_element(const std::vector<float>& vec);

/**
 * @brief Selects n random elements from a vector.
 * @param vec Input vector.
 * @param n Number of elements to select.
 * @return Vector of selected elements.
 */
std::vector<int> select_random_elements(const std::vector<int>& vec, int n);

/**
 * @brief Computes the cumulative distribution function (CDF) of a discrete distribution.
 * @param probabilities Input probability distribution.
 * @return CDF of the input distribution.
 */
std::vector<float> computeCDF(const std::vector<float>& probabilities);

/**
 * @brief Generates a random value from a discrete distribution using its CDF.
 * @param probabilities Input probability distribution.
 * @return Randomly generated value.
 */
int generateRandomDiscrete(const std::vector<float>& probabilities);

/**
 * @brief Normalizes a distribution into a probability distribution.
 * @param values Input values.
 * @return Normalized probability distribution.
 */
std::vector<float> normalizeDistribution(const std::vector<float>& values);
