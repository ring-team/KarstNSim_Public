/**
 * @file    simplex_noise.h
 * @brief   A Perlin Simplex Noise C++ Implementation (1D, 2D, 3D).
 *
 * Copyright (c) 2014-2018 Sebastien Rombauts (sebastien.rombauts@gmail.com)
 *
 * Distributed under the MIT License (MIT) (See accompanying file LICENSE.txt
 * or copy at http://opensource.org/licenses/MIT)
 */
#pragma once

#include <cstddef>  // size_t
#include <cstdint>
#include <algorithm>
#include <array>
#include <random>
 /**
  * @brief A Perlin Simplex Noise C++ Implementation (1D, 2D, 3D, 4D).
  */
class SimplexNoise {
public:
	//! Permutation table used to hash lattice coordinates. Each instance owns its own table.
	using Permutation = std::array<std::uint8_t, 256>;

	// 1D Perlin simplex noise
	float noise(float x) const;
	// 2D Perlin simplex noise
	float noise(float x, float y) const;
	// 3D Perlin simplex noise
	float noise(float x, float y, float z) const;

	// Fractal/Fractional Brownian Motion (fBm) noise summation
	float fractal(size_t octaves, float x) const;
	float fractal(size_t octaves, float x, float y) const;
	float fractal(size_t octaves, float x, float y, float z) const;

	//initialize seed coherent permutation table (shuffles the reference table with a copy of rng)
	void initialize_permutation_table(std::mt19937 rng);
	//replace this instance's permutation table
	void set_permutation(const Permutation& permutation) { mPerm = permutation; }
	//permutation table of this instance
	const Permutation& permutation() const { return mPerm; }
	//reference (unshuffled) Perlin permutation table
	static const Permutation& reference_permutation();
	int print_perm();

	/**
	 * Constructor of to initialize a fractal noise summation
	 *
	 * @param[in] frequency    Frequency ("width") of the first octave of noise (default to 1.0)
	 * @param[in] amplitude    Amplitude ("height") of the first octave of noise (default to 1.0)
	 * @param[in] lacunarity   Lacunarity specifies the frequency multiplier between successive octaves (default to 2.0).
	 * @param[in] persistence  Persistence is the loss of amplitude between successive octaves (usually 1/lacunarity)
	 */
	explicit SimplexNoise(float frequency = 1.0f,
		float amplitude = 1.0f,
		float lacunarity = 2.0f,
		float persistence = 0.5f) :
		mFrequency(frequency),
		mAmplitude(amplitude),
		mLacunarity(lacunarity),
		mPersistence(persistence),
		mPerm(reference_permutation()) {


	}

private:
	// Hash an integer using this instance's permutation table
	std::uint8_t hash(std::int32_t i) const {
		return mPerm[static_cast<std::uint8_t>(i)];
	}

	// Parameters of Fractional Brownian Motion (fBm) : sum of N "octaves" of noise
	float mFrequency;   ///< Frequency ("width") of the first octave of noise (default to 1.0)
	float mAmplitude;   ///< Amplitude ("height") of the first octave of noise (default to 1.0)
	float mLacunarity;  ///< Lacunarity specifies the frequency multiplier between successive octaves (default to 2.0).
	float mPersistence; ///< Persistence is the loss of amplitude between successive octaves (usually 1/lacunarity)
	Permutation mPerm; ///< Permutation table owned by this instance
};
