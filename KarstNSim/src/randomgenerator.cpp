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

double generateNormalRandom(double mean, double stddev) {
	std::normal_distribution<double> normalDistribution(mean, stddev);
	return normalDistribution(globalRng);
}

int generateRandomInt(int min, int max) {
	std::uniform_int_distribution<int> distribution(min, max);
	return distribution(globalRng);
}

double generateRandomDouble(double min, double max) {
	std::uniform_real_distribution<double> distribution(min, max);
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

std::vector<double> generateGaussianVector(std::size_t size, double mean, double stddev) {
	std::vector<double> result;
	result.reserve(size);

	std::normal_distribution<double> distribution(mean, stddev);
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

std::vector<double> generateRandomDoubleVector(std::size_t size, double minValue, double maxValue) {
	std::vector<double> result;
	result.reserve(size);

	std::uniform_real_distribution<double> distribution(minValue, maxValue);
	for (std::size_t i = 0; i < size; ++i) {
		result.push_back(distribution(globalRng));
	}

	return result;
}
