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

// Declare the random number generator and seed globally
extern std::mt19937 globalRng;
extern std::vector<std::uint32_t> globalSeed;

// Function to initialize the RNG with a seed
void initializeRng(const std::vector<std::uint32_t>& seed);

// Function to generate a random number from a normal distribution
double generateNormalRandom(double mean, double stddev);

// Function to generate a random integer within a range
int generateRandomInt(int min, int max);

// Function to generate a random double within a range
double generateRandomDouble(double min, double max);

// Function to generate a random boolean
bool generateRandomBool();

// Function to generate a random character
char generateRandomChar();

// Function to generate a random index (useful for picking a random element from a container)
std::size_t generateRandomIndex(std::size_t containerSize);

// Function to shuffle a container using globalRng
template <typename ForwardIt>
void shuffleContainer(ForwardIt first, ForwardIt last) {
	std::shuffle(first, last, globalRng);
}

// Function to generate a random Gaussian (normal) distribution vector
std::vector<double> generateGaussianVector(std::size_t size, double mean, double stddev);

// Function to generate a random vector of integers
std::vector<int> generateRandomIntVector(std::size_t size, int minValue, int maxValue);

// Function to generate a random vector of doubles
std::vector<double> generateRandomDoubleVector(std::size_t size, double minValue, double maxValue);
