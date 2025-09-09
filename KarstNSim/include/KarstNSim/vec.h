/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr for some generic functions and modifications to other classes
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

Copyright (c) 2021 Axel Paris
Author : Axel Paris, for classes Vector2i, Vector2, Vector3
If you use this code, pleace cite : Paris et al., 2021, Computer Graphic Forum.

***************************************************************/

#pragma once

/**
@file vec.h
@brief Utility functions, and classes for simple vector objects
@authors Axel PARIS, Augustin GOUY
**/

#include <KarstNSim/nanoflann.hpp>
#include <cmath>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <iterator>

/* Forward Declarations */

/*!
\struct Vector2i
\brief A structure representing a 2D vector with integer components.
This structure is typically used to store coordinates or vectors in a 2D space, where the components are integers.
*/
struct Vector2i;

/*!
\struct Vector2
\brief A structure representing a 2D vector with floating-point components.
This structure is typically used to store coordinates or vectors in a 2D space, where the components are floating-point numbers for increased precision.
*/
struct Vector2;

/*!
\struct Vector3
\brief A structure representing a 3D vector with floating-point components.
This structure is typically used to store coordinates or vectors in a 3D space, with floating-point components for precision.
*/
struct Vector3;

// Maths utility
namespace KarstNSim
{
	/*!
	\brief Clamps a value between two limits
	\param x Value
	\param a, b Interval values.
	\return Real.
	*/
	inline float Clamp(float x, float a = 0.0f, float b = 1.0f)
	{
		return x < a ? a : x > b ? b : x;
	}

	/*!
	\brief Create a linear step.
	\param x Value
	\param a, b Interval values.
	\return Real in unit inverval.
	*/
	inline float Step(const float& x, const float& a, const float& b)
	{
		if (x < a)
		{
			return 0.0f;
		}
		else if (x > b)
		{
			return 1.0f;
		}
		else
		{
			return (x - a) / (b - a);
		}
	}

	template<typename T>
	inline T Min(T a, T b)
	{
		return a < b ? a : b;
	}

	template<typename T>
	inline T Max(T a, T b)
	{
		return a > b ? a : b;
	}

	template<typename T>
	inline T Lerp(T a, T b, float t)
	{
		return (a * (1.0f - t)) + (b * t);
	}

	inline float Abs(float a)
	{
		return a < 0 ? -a : a;
	}

	inline float CubicSmoothCompact(float x, float r)
	{
		return (x > r) ? 0.0f : (1.0f - x / r) * (1.0f - x / r) * (1.0f - x / r);
	}
	template <class T> T _square(T a) { return a * a; };

	inline float CubicSigmoid(float x, float r, float t)
	{
		if (x > 0.0f)
		{
			if (x < r)
			{
				return x * (1.0f + x * ((3.0f * t - 2.0f * r) / (r * r) + x * (r - 2.0f * t) / (r * r * r)));
			}
			else
			{
				return t;
			}
		}
		else
		{
			if (x > -r)
			{
				// Use symmetric
				float y = -x;
				return  -(y * (1.0f + y * ((3.0f * t - 2.0f * r) / (r * r) + y * (r - 2.0f * t) / (r * r * r))));
			}
			else
			{
				return -t;
			}
		}
	}

	inline float CubicSmooth(float x)
	{
		return x * x * (3.0 - 2.0 * x);
	}

	inline float CubicSmooth(float x, float r)
	{
		return (1.0 - x / r) * (1.0 - x / r) * (1.0 - x / r);
	}

	inline float CubicSmoothStep(float x, float a, float b)
	{
		if (x < a)
			return 0.0;
		else if (x > b)
			return 1.0;
		else
			return 1.0 - CubicSmooth((x - a) * (x - a), (b - a) * (b - a));
	}

	inline float QuinticSmooth(float t)
	{
		return pow(t, 3.0) * (t * (t * 6.0 - 15.0) + 10.0);
	}

	inline float Pow(float x, float e)
	{
		if (x == 0.0)
		{
			return 0.0;
		}
		else if (x > 0.0)
		{
			return pow(x, e);
		}
		else
		{
			return -pow(-x, e);
		}
	}

	inline float calculateCircumradius(float a, float b, float c) {
		float s = (a + b + c) / 2;
		float circumradius = (a * b * c) / (4 * sqrt(s * (s - a) * (s - b) * (s - c)));
		return circumradius;
	}

	/*!
	\brief Checks whether a triangle with given side lengths is valid.

	This function verifies the validity of a triangle based on the given side lengths, ensuring the triangle inequality holds and the triangle has a non-zero area within a relative precision.

	\param a Length of the first side of the triangle.
	\param b Length of the second side of the triangle.
	\param c Length of the third side of the triangle.
	\return True if the triangle is valid; false otherwise.

	### Conditions for validity:
	1. The side lengths must be non-negative.
	2. The triangle inequality must hold: \( a + b > c \), \( a + c > b \), \( b + c > a \), considering a relative epsilon to handle floating-point precision.
	3. The area of the triangle, calculated using Heron's formula, must be greater than a small epsilon to avoid degenerate cases.

	### Notes:
	- The function uses a relative epsilon based on the magnitudes of the side lengths to ensure robustness against floating-point precision issues.
	- A degenerate triangle (with zero or nearly zero area) will be considered invalid.
	*/
	inline bool isTriangleValid(float a, float b, float c) {
		// Choose a small multiplier for epsilon (e.g., 1e-5)
		const float epsilonMultiplier = 1e-5f;

		// Calculate a relative epsilon based on the magnitudes of the side lengths
		float epsilon = epsilonMultiplier * std::max({ std::abs(a), std::abs(b), std::abs(c) });

		// Check for non-negative side lengths
		if (a < 0 || b < 0 || c < 0) {
			return false;
		}

		// Check the Triangle Inequality with the relative epsilon
		if (a + b > c + epsilon && a + c > b + epsilon && b + c > a + epsilon) {
			// Check for degenerate cases (zero area) with the relative epsilon
			float s = (a + b + c) / 2;
			float area = sqrt(s * (s - a) * (s - b) * (s - c));
			return area > epsilon;
		}

		return false;
	}
}


/* Vector2i */
struct Vector2i
{
public:
	int x, y;

	explicit Vector2i() : x(0), y(0) {}
	explicit Vector2i(int n) : x(n), y(n) {}
	explicit Vector2i(int x, int y) : x(x), y(y) {}

	Vector2i operator-(const Vector2i& u) const
	{
		return Vector2i(x - u.x, y - u.y);
	}
	Vector2i operator+(const Vector2i& u) const
	{
		return Vector2i(x + u.x, y + u.y);
	}
};


/* Vector3 */
struct Vector3
{
public:
	float x, y, z;

	Vector3() : x(0.0f), y(0.0f), z(0.0f) { }
	Vector3(float n) : x(n), y(n), z(n) {}
	Vector3(float x, float y, float z) : x(x), y(y), z(z) {}
	Vector3(const Vector3 &v) = default;
	Vector3& operator= (const Vector3&) = default;

	friend bool operator> (const Vector3&, const Vector3&);
	friend bool operator< (const Vector3&, const Vector3&);
	friend bool operator>= (const Vector3&, const Vector3&);
	friend bool operator<= (const Vector3&, const Vector3&);



	Vector3& operator+= (const Vector3& v)
	{
		x += v.x;
		y += v.y;
		z += v.z;
		return *this;
	}

	//Vector3& operator=(const Vector3 &v)
	//{
	//	return *this;
	//}

	Vector3 operator-= (const Vector3& v)
	{
		return Vector3(x - v.x, y - v.y, z - v.z);
	}
	Vector3 operator*= (float f)
	{
		return Vector3(x * f, y * f, z * f);
	}
	Vector3 operator/= (float f)
	{
		return Vector3(x / f, y / f, z / f);
	}
	Vector3 operator*(const Vector3& u) const
	{
		return Vector3(x * u.x, y * u.y, z * u.z);
	}
	Vector3 operator*(float k) const
	{
		return Vector3(x * k, y * k, z * k);
	}
	Vector3 operator/(float k) const
	{
		return Vector3(x / k, y / k, z / k);
	}
	bool operator==(const Vector3& u) const
	{
		float eps = 1e-2f;
		return (std::abs(x - u.x) < eps && std::abs(y - u.y) < eps && std::abs(z - u.z) < eps);
	}
	bool operator!=(const Vector3& u) const
	{
		float eps = 1e-2f;
		return (std::abs(x - u.x) > eps || std::abs(y - u.y) > eps || std::abs(z - u.z) > eps);
	}
	Vector3 operator-(const Vector3& u) const
	{
		return Vector3(x - u.x, y - u.y, z - u.z);
	}
	Vector3 operator+(const Vector3& u) const
	{
		return Vector3(x + u.x, y + u.y, z + u.z);
	}
	Vector3 operator+(float k) const
	{
		return Vector3(x + k, y + k, z + k);
	}
	float operator[](int i) const
	{
		if (i == 0)
			return x;
		else if (i == 1)
			return y;
		return z;
	}
	float& operator[](int i)
	{
		if (i == 0)
			return x;
		else if (i == 1)
			return y;
		return z;
	}

	friend std::ostream& operator<<(std::ostream& stream, const Vector3& u);
	inline float Max() const
	{
		return KarstNSim::Max(KarstNSim::Max(x, y), z);
	}
	inline float Min() const
	{
		return KarstNSim::Min(KarstNSim::Min(x, y), z);
	}
	inline int MaxIndex() const
	{
		if (x >= y)
		{
			if (x >= z)
				return 0;
			else
				return 2;
		}
		else
		{
			if (y >= z)
				return 1;
			else
				return 2;
		}
	}
	static inline Vector3 Min(const Vector3& a, const Vector3& b)
	{
		return Vector3(KarstNSim::Min(a.x, b.x), KarstNSim::Min(a.y, b.y), KarstNSim::Min(a.z, b.z));
	}
	static inline Vector3 Max(const Vector3& a, const Vector3& b)
	{
		return Vector3(KarstNSim::Max(a.x, b.x), KarstNSim::Max(a.y, b.y), KarstNSim::Max(a.z, b.z));
	}
};

namespace KarstNSim
{
	inline std::ostream& operator<<(std::ostream& stream, const Vector3& u)
	{
		stream << "(" << u.x << ", " << u.y << ", " << u.z << ");";
		return stream;
	}
	inline Vector3 Cross(const Vector3& u, const Vector3& v)
	{
		return Vector3((u.y * v.z) - (u.z * v.y), (u.z * v.x) - (u.x * v.z), (u.x * v.y) - (u.y * v.x));
	}

	inline float Dot2D(const Vector3& a, const Vector3& b) {
		return a.x * b.x + a.y * b.y;
	}

	inline float Dot(const Vector3& u, const Vector3& v)
	{
		return u.x * v.x + u.y * v.y + u.z * v.z;
	}

	// reorder vector with given vector of indices
	inline void reorder(std::vector<int> &vA, std::vector<int> vOrder)
	{
		assert(vA.size() == vOrder.size());
		// for all elements to put in place
		for (int i = 0; i < vA.size() - 1; ++i)
		{
			// while the element i is not yet in place 
			while (i != vOrder[i])
			{
				// swap it with the element at its final place
				int alt = vOrder[i];
				std::swap(vA[i], vA[alt]);
				std::swap(vOrder[i], vOrder[alt]);
			}
		}
	}

	inline float magnitude2D(const Vector3& u)
	{
		return sqrt(u.x * u.x + u.y * u.y);
	}

	inline float squaredmagnitude2D(const Vector3& u)
	{
		return u.x * u.x + u.y * u.y;
	}

	inline float magnitude(const Vector3& u)
	{
		return sqrt(u.x * u.x + u.y * u.y + u.z * u.z);
	}

	inline float squaredmagnitude(const Vector3& u)
	{
		return u.x * u.x + u.y * u.y + u.z * u.z;
	}

	inline Vector3 mid(const Vector3& u, const Vector3& v)
	{
		return Vector3((u.x + v.x) / 2, (u.y + v.y) / 2, (u.z + v.z) / 2);
	}

	inline void squaredistance_to_segment(const Vector3& pt, const Vector3& p1, const Vector3& p2, float& squaredDistance, Vector3& closestPoint)
	{
		Vector3 v = p2 - p1;
		Vector3 w = pt - p1;
		float eps = 1e-5f;

		float c1 = Dot(w, v);
		if (c1 <= eps) {
			squaredDistance = squaredmagnitude(w);
			closestPoint = p1;
			return;
		}

		float c2 = Dot(v, v);
		if (c2 <= c1) {
			squaredDistance = squaredmagnitude(pt - p2);
			closestPoint = p2;
			return;
		}

		float b = c1 / c2;
		Vector3 pb = p1 + v * b;
		squaredDistance = squaredmagnitude(pt - pb);
		closestPoint = pb;
	}

	inline void squaredistance_to_segment2D(const Vector3& pt, const Vector3& p1, const Vector3& p2, float& squaredDistance, Vector3& closestPoint)
	{
		Vector3 v = p2 - p1;
		Vector3 w = pt - p1;
		float eps = 1e-5f;

		float c1 = Dot2D(w, v);
		if (c1 <= eps) {
			squaredDistance = squaredmagnitude2D(w);
			closestPoint = p1;
			return;
		}

		float c2 = Dot2D(v, v);
		if (c2 <= c1) {
			squaredDistance = squaredmagnitude2D(pt - p2);
			closestPoint = p2;
			return;
		}

		float b = c1 / c2;
		Vector3 pb = p1 + v * b;
		squaredDistance = squaredmagnitude2D(pt - pb);
		closestPoint = pb;
	}


	inline float magnitudeNormalized(const Vector3& u, float du, float dv, float dw)
	{
		return sqrt((u.x * u.x) / (du*du) + (u.y * u.y) / (dv*dv) + (u.z * u.z) / (dw*dw));
	}

	inline float squaredmagnitudeNormalized(const Vector3& u, float one_over_delta_x_squared, float one_over_delta_y_squared, float one_over_delta_z_squared)
	{
		return (u.x * u.x)*one_over_delta_x_squared + (u.y * u.y)*one_over_delta_y_squared + (u.z * u.z)*one_over_delta_z_squared;
	}

	inline Vector3 Normalize(const Vector3& v)
	{
		float magSquared = v.x * v.x + v.y * v.y + v.z * v.z;

		// Avoid dividing by zero by checking if magnitude squared is non-zero
		return (magSquared > 0.0) ? v * (1.0 / sqrt(magSquared)) : Vector3(0.0, 0.0, 0.0);
	}
	inline Vector3 operator-(const Vector3& v)
	{
		return Vector3(-v.x, -v.y, -v.z);
	}
	inline bool operator>(const Vector3& u, const Vector3& v)
	{
		return (u.x > v.x) && (u.y > v.y) && (u.z > v.z);
	}
	inline bool operator<(const Vector3& u, const Vector3& v)
	{
		return (u.x < v.x) && (u.y < v.y) && (u.z < v.z);
	}
	inline bool operator>=(const Vector3& u, const Vector3& v)
	{
		return (u.x >= v.x) && (u.y >= v.y) && (u.z >= v.z);
	}
	inline bool operator<=(const Vector3& u, const Vector3& v)
	{
		return (u.x <= v.x) && (u.y <= v.y) && (u.z <= v.z);
	}
	inline Vector3 operator*(float a, const Vector3& v)
	{
		return v * a;
	}
	inline Vector3 Abs(const Vector3& u)
	{
		return Vector3(u[0] > 0.0 ? u[0] : -u[0], u[1] > 0.0 ? u[1] : -u[1], u[2] > 0.0 ? u[2] : -u[2]);
	}

	/*!
	\brief Standardizes a vector of values to a specified range.

	This function rescales the values in the vector to fit within the range [min_value, max_value].
	No-data values (normally -99999) are excluded from the standardization process.

	\param vec Vector of values to be standardized.
	\param min_value Minimum value of the target range.
	\param max_value Maximum value of the target range.

	\details
	- The function finds the minimum and maximum values in the vector (ignoring no-data values).
	- It applies min-max scaling to the valid values.
	- If all values are the same, they are set to the midpoint of the target range.
	*/
	inline void standardize_to_range(std::vector<float>& vec, float min_value, float max_value) {

		// Find the minimum and maximum values in the vector, ignoring the no-data values
		float min_val = std::numeric_limits<float>::max();
		float max_val = std::numeric_limits<float>::lowest();

		for (const auto& val : vec) {
			if (val >= -10000) { // Exclude no-data values (= normally -99999) based on epsilon
				min_val = std::min(min_val, val);
				max_val = std::max(max_val, val);
			}
		}

		// Compute the range of the vector values
		float range = max_val - min_val;

		// Iterate through each element in the vector
		for (int i = 0; i < vec.size(); ++i) {
			// Exclude no-data values from standardization
			if (vec[i] >= -10000) { // -99999
				// Apply min-max scaling to standardize the value between 0 and 1
				if (range != 0) {
					vec[i] = ((vec[i] - min_val) / range) * (max_value - min_value) + min_value;
				}
				else {
					// If range is 0, set all elements to the mid-point of the range
					vec[i] = (max_value + min_value) / 2.0;
				}
			}
		}
	}


	/*!
	\brief Converts a 3D point in a grid to a flattened 1D index.

	This operation is known as raveling, mapping a 3D point to a 1D array index.

	\param dim1 Number of divisions in the x-dimension.
	\param dim2 Number of divisions in the y-dimension.
	\param dim3 Number of divisions in the z-dimension.
	\param point 3D point with normalized (0,1) coordinates.
	\return Flattened index corresponding to the input point.
	*/
	inline int64_t flatten_indices(int dim1, int dim2, int dim3, const Vector3& point) {

		int64_t flattened_index = 0;

		int64_t x = static_cast<int64_t>(std::floor(point.x * dim1));
		int64_t y = static_cast<int64_t>(std::floor(point.y * dim2));
		int64_t z = static_cast<int64_t>(std::floor(point.z * dim3));

		// Ensure that the indices are within bounds
		x = std::max<int64_t>(0, std::min<int64_t>(x, dim1 - 1));
		y = std::max<int64_t>(0, std::min<int64_t>(y, dim2 - 1));
		z = std::max<int64_t>(0, std::min<int64_t>(z, dim3 - 1));

		flattened_index = x + y * dim1 + z * dim1 * dim2;
		return flattened_index;
	}


	/*!
	\brief Converts a 3D point to a flattened index for a grid with custom scaling.

	Similar to flatten_indices but allows for grids with arbitrary dimensions and scaling factors.

	\param dim1 Number of divisions in the x-dimension.
	\param dim2 Number of divisions in the y-dimension.
	\param dim3 Number of divisions in the z-dimension.
	\param dimx Physical size of the grid in the x-dimension.
	\param dimy Physical size of the grid in the y-dimension.
	\param dimz Physical size of the grid in the z-dimension.
	\param point 3D point with normalized (0,1) coordinates.
	\return Flattened index corresponding to the input point.
	*/
	inline int64_t flatten_indices2(const int& dim1, const int& dim2, const int& dim3, const float& dimx, const float& dimy, const float& dimz, const Vector3& point) {

		int64_t flattened_index = 0;

		int64_t x = static_cast<int64_t>(std::floor(point.x * dim1 / dimx));
		int64_t y = static_cast<int64_t>(std::floor(point.y * dim2 / dimy));
		int64_t z = static_cast<int64_t>(std::floor(point.z * dim3 / dimz));

		// Ensure that the indices are within bounds
		x = std::max<int64_t>(0, std::min<int64_t>(x, dim1 - 1));
		y = std::max<int64_t>(0, std::min<int64_t>(y, dim2 - 1));
		z = std::max<int64_t>(0, std::min<int64_t>(z, dim3 - 1));

		flattened_index = x + y * dim1 + z * static_cast<int64_t>(dim1) * dim2;

		return flattened_index;
	}

	/*!
	\brief Converts a flattened 1D index back to a 3D coordinate.

	This operation is known as unraveling, mapping a 1D array index to 3D grid coordinates.

	\param idx_in Flattened index in a 1D array.
	\param intdim1 Number of divisions in the x-dimension.
	\param intdim2 Number of divisions in the y-dimension.
	\param index_init_norm Determines if indices start at 0 (default) or 1.
	\return 3D coordinates corresponding to the input index.
	*/
	inline Vector3 unravel3D(const int64_t& idx_in, const int& intdim1, const int& intdim2, int index_init_norm = 0) {
		int64_t i, j, k;

		if (index_init_norm == 0) {
			k = idx_in / static_cast<int64_t>(intdim1 * intdim2);
			j = (idx_in - k * static_cast<int64_t>(intdim1) * intdim2) / intdim1;
			i = idx_in - k * static_cast<int64_t>(intdim1) * intdim2 - j * intdim1;
		}
		else {
			k = idx_in / static_cast<int64_t>(intdim1 * intdim2) + 1;
			j = (idx_in - (k - 1) * static_cast<int64_t>(intdim1) * intdim2) / intdim1 + 1;
			i = idx_in - (k - 1) * static_cast<int64_t>(intdim1) * intdim2 - (j - 1) * intdim1;
		}

		return Vector3(i, j, k);
	}

	/*!
	\brief Converts a 1D index from one grid to another, accounting for differing grid sizes.

	\param dim1_in Number of divisions in the x-dimension of the input grid.
	\param dim2_in Number of divisions in the y-dimension of the input grid.
	\param dim3_in Number of divisions in the z-dimension of the input grid.
	\param dim1_out Number of divisions in the x-dimension of the output grid.
	\param dim2_out Number of divisions in the y-dimension of the output grid.
	\param dim3_out Number of divisions in the z-dimension of the output grid.
	\param idx_in Flattened index in the input grid.
	\param index_init_norm Determines if indices start at 0 (default) or 1.
	\return Flattened index in the output grid.
	*/
	inline int64_t convert_index(int dim1_in, int dim2_in, int dim3_in, int dim1_out, int dim2_out, int dim3_out, int64_t idx_in, int index_init_norm = 0) {
		Vector3 coord = unravel3D(idx_in, dim1_in, dim2_in, index_init_norm);
		int64_t i_out = static_cast<int64_t>(std::floor(coord.x * dim1_out / dim1_in));
		int64_t j_out = static_cast<int64_t>(std::floor(coord.y * dim2_out / dim2_in));
		int64_t k_out = static_cast<int64_t>(std::floor(coord.z * dim3_out / dim3_in));
		int64_t idx_out = i_out + j_out * dim1_out + k_out * dim1_out * dim2_out;
		return idx_out;
	}

	/*!
	\brief Scans a local region around a 3D point within a grid for nearby elements.

	The method identifies grid elements within a specified radius from a given point.

	\param r Search radius.
	\param r_min Minimum resolution of the grid.
	\param point 3D point in normalized coordinates (0,1).
	\param dim Number of divisions in each dimension of the grid.
	\return List of flattened indices representing nearby elements.
	*/
	inline std::vector<int64_t> scan_area(float r, float r_min, const Vector3& point, int dim) {
		int64_t flattened_index = flatten_indices(dim, dim, dim, point);
		Vector3 coord = unravel3D(flattened_index, dim, dim);
		std::vector<int64_t> scanned_indices;
		int region_size = static_cast<int>(std::ceil(r * std::sqrt(3) / r_min)); // we're looking to scan all the elements at a distance less than region_size of the index
		int64_t x_i, y_i, z_i;
		for (x_i = -region_size; (x_i <= region_size); x_i++) {
			if ((0 <= coord.x + x_i) && (coord.x + x_i < dim)) { // we also check that those elements aren't outside the grid
				for (y_i = -region_size; (y_i <= region_size); y_i++) {
					if ((0 <= coord.y + y_i) && (coord.y + y_i < dim)) {
						for (z_i = -region_size; (z_i <= region_size); z_i++) {
							if ((0 <= coord.z + z_i) && (coord.z + z_i < dim)) {
								scanned_indices.push_back(flattened_index + x_i + y_i * dim + z_i * dim * dim);
							}
						}
					}
				}
			}
		}
		return scanned_indices;
	}

	/*!
	\brief Scans a local region around a 3D point in a grid with arbitrary dimensions.

	Similar to scan_area but supports grids with arbitrary scaling and dimensions.

	\param r Search radius.
	\param r_min Minimum resolution of the grid.
	\param point 3D point in normalized coordinates (0,1).
	\param dim1 Number of divisions in the x-dimension.
	\param dim2 Number of divisions in the y-dimension.
	\param dim3 Number of divisions in the z-dimension.
	\param dimx Physical size of the grid in the x-dimension.
	\param dimy Physical size of the grid in the y-dimension.
	\param dimz Physical size of the grid in the z-dimension.
	\return List of flattened indices representing nearby elements.
	*/
	inline std::vector<int64_t> scan_area2(const float& r, const float& r_min, const Vector3& point, const int& dim1, const int& dim2,
		const int& dim3, const float& dimx, const float& dimy, const float& dimz) {

		int64_t flattened_index = flatten_indices2(dim1, dim2, dim3, dimx, dimy, dimz, point);
		Vector3 coord = unravel3D(flattened_index, dim1, dim2);
		std::vector<int64_t> scanned_indices;
		int region_size = static_cast<int>(std::ceil(r * std::sqrt(3) / r_min)); // we're looking to scan all the elements at a distance less than region_size of the index
		int64_t x_i, y_i, z_i;

		for (x_i = -region_size; (x_i <= region_size); x_i++) {
			if ((0 <= coord.x + x_i) && (coord.x + x_i < dim1)) { // we also check that those elements aren't outside the grid
				for (y_i = -region_size; (y_i <= region_size); y_i++) {
					if ((0 <= coord.y + y_i) && (coord.y + y_i < dim2)) {
						for (z_i = -region_size; (z_i <= region_size); z_i++) {
							if ((0 <= coord.z + z_i) && (coord.z + z_i < dim3)) {
								scanned_indices.push_back(flattened_index + x_i + y_i * dim1 + z_i * dim1 * dim2);
							}
						}
					}
				}
			}
		}
		return scanned_indices;
	}

	/*!
	\brief Checks if a 2D point lies inside a triangle IN MAP VIEW (2D, so no z check!!)

	\param A First vertex of the triangle.
	\param B Second vertex of the triangle.
	\param C Third vertex of the triangle.
	\param P Point to check.
	\param tolerance Tolerance for inclusion, default is 1e-2.
	\return True if the point is inside the triangle, false otherwise.
	*/
	inline bool isInsideTriangle(const Vector3& A, const Vector3& B, const Vector3& C, const Vector3& P, const float tolerance = 1e-2) {
		// Calculate barycentric coordinates
		const float detT = (B.y - C.y) * (A.x - C.x) + (C.x - B.x) * (A.y - C.y);
		const float alpha = ((B.y - C.y) * (P.x - C.x) + (C.x - B.x) * (P.y - C.y)) / detT;
		const float beta = ((C.y - A.y) * (P.x - C.x) + (A.x - C.x) * (P.y - C.y)) / detT;
		const float gamma = 1.0 - alpha - beta;

		// Check if the point is inside the triangle with tolerance
		return alpha >= -tolerance && beta >= -tolerance && gamma >= -tolerance;
	}

	/*!
	\brief Adds a Vector3 to a list if it is not already present within a tolerance.

	\param vec Vector of Vector3 objects to check.
	\param newVector The new Vector3 to add.
	\param tolerance Tolerance for determining equality, default is 1e-5.
	*/
	inline void pushBackIfNotPresent(std::vector<Vector3>& vec, const Vector3& newVector, float tolerance = 1e-5) {
		if (std::find_if(vec.begin(), vec.end(), [&](const Vector3& v) {
			return std::abs(v.x - newVector.x) < tolerance && std::abs(v.y - newVector.y) < tolerance &&
				std::abs(v.z - newVector.z) < tolerance;
		}) == vec.end()) {
			vec.push_back(newVector);
		}
	}

	inline void vect_product(Vector3 vec1, Vector3 vec2, Vector3& result) {
		result.x = vec1.y*vec2.z - vec1.z*vec2.y;
		result.y = vec1.z*vec2.x - vec1.x*vec2.z;
		result.z = vec1.x*vec2.y - vec1.y*vec2.x;
	}

	inline float scalar_product(const Vector3& v1, const Vector3& v2) {
		return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
	}

	/*!
	\brief Checks if a point is below a triangle in 3D space.

	This function determines whether a point lies below a triangle defined by three vertices
	by checking its signed distance to the plane of the triangle.

	\param new_horizon_pt Point to test.
	\param P0 First vertex of the triangle.
	\param P1 Second vertex of the triangle.
	\param P2 Third vertex of the triangle.
	\param tolerance Tolerance for determining coplanarity, default is 1.
	\return True if the point is below the triangle or coplanar within the given tolerance, false otherwise.

	\details
	- The function computes the normal vector of the triangle's plane and its signed distance to the point.
	- A point is considered below the triangle if:
	  - The triangle's normal points towards the positive z-axis and the signed distance is less than `-tolerance`.
	  - The triangle's normal points towards the negative z-axis and the signed distance is greater than `tolerance`.
	*/
	inline bool is_point_under_triangle(const Vector3& new_horizon_pt, const Vector3& P0, const Vector3& P1, const Vector3& P2, const float tolerance = 1) {
		Vector3 edge1 = { P1.x - P0.x, P1.y - P0.y, P1.z - P0.z };
		Vector3 edge2 = { P2.x - P0.x, P2.y - P0.y, P2.z - P0.z };
		Vector3 normal;
		vect_product(edge1, edge2, normal);
		float d = -scalar_product(normal, P0);
		float res = scalar_product(new_horizon_pt, normal) + d;

		// Check if the triangle is facing towards positive z-axis
		bool isFacingPositiveZ = normal.z > 0;

		if (std::abs(res) <= tolerance) {
			return true;  // Point is coplanar or very close
		}
		else if ((isFacingPositiveZ && res < -tolerance) || (!isFacingPositiveZ && res > tolerance)) {
			return true;  // Point is below the triangle
		}
		else {
			return false;  // Point is above the triangle
		}
	}


	/*!
	\brief Computes the orthogonal projection distance from a point to the plane of a triangle.

	This function calculates the absolute value of the signed distance from a point to the plane
	defined by three vertices of a triangle.

	\param pt Point to test.
	\param A First vertex of the triangle.
	\param B Second vertex of the triangle.
	\param C Third vertex of the triangle.
	\return Absolute value of the orthogonal projection distance from the point to the plane.

	\details
	- The function computes the normal vector of the triangle's plane.
	- The signed distance is calculated as the dot product of the normalized normal vector and the vector from one triangle vertex to the point.
	*/
	inline float barycentric(const Vector3& pt, const Vector3& A, const Vector3& B, const Vector3& C) {

		// Compute the normal vector of the plane defined by the triangle
		Vector3 normal = Cross(B - A, C - A);

		// Normalize the normal vector
		const float normalLength = std::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
		normal.x /= normalLength;
		normal.y /= normalLength;
		normal.z /= normalLength;

		// Compute the signed distance to the plane
		const float signedDistance = Dot(normal, pt - A);

		// Return the absolute value of the signed distance (orthogonal projection distance)
		return std::abs(signedDistance);
	}



	/*!
	\brief Computes the signed distance from a point to the plane defined by a triangle.

	This function calculates the signed distance between a point and the plane determined by three triangle vertices.
	The sign of the distance indicates whether the point lies above or below the plane.

	\param pt Point to test.
	\param A First vertex of the triangle.
	\param B Second vertex of the triangle.
	\param C Third vertex of the triangle.
	\return Signed distance from the point to the plane.

	\details
	- The function computes the normal vector of the triangle's plane by taking the cross product of two edges.
	- The signed distance is calculated as the dot product of the normalized normal vector and the vector from one triangle vertex to the point.
	*/
	inline float signed_distance_to_plane(const Vector3& pt, const Vector3& A, const Vector3& B, const Vector3& C) {
		const Vector3 normal = Normalize(Cross((B - A), (C - A)));
		return Dot(normal, pt - A);
	}
}

/* Vector2 */
struct Vector2
{
public:
	float x, y;

	explicit Vector2() : x(0.0f), y(0.0f) { }
	explicit Vector2(const Vector2i& v) : x(float(v.x)), y(float(v.y)) { }
	explicit Vector2(float n) : x(n), y(n) { }
	explicit Vector2(float x, float y) : x(x), y(y) { }
	explicit Vector2(const Vector3& v) : x(v.x), y(v.z) { }

	friend bool operator> (const Vector2&, const Vector2&);
	friend bool operator< (const Vector2&, const Vector2&);
	friend bool operator>= (const Vector2&, const Vector2&);
	friend bool operator<= (const Vector2&, const Vector2&);

	Vector2 operator+= (const Vector2& v)
	{
		return Vector2(x + v.x, y + v.y);
	}
	Vector2 operator-= (const Vector2& v)
	{
		return Vector2(x - v.x, y - v.y);
	}
	Vector2 operator*= (float f)
	{
		return Vector2(x * f, y * f);
	}
	Vector2 operator/= (float f)
	{
		return Vector2(x / f, y / f);
	}
	Vector2 operator*(const Vector2& v) const
	{
		return Vector2(x * v.x, y * v.y);
	}
	Vector2 operator*(float k) const
	{
		return Vector2(x * k, y * k);
	}
	Vector2 operator/(float k) const
	{
		return Vector2(x / k, y / k);
	}
	bool operator==(const Vector2& u) const
	{
		return (x == u.x && y == u.y);
	}
	Vector2 operator-(const Vector2& u) const
	{
		return Vector2(x - u.x, y - u.y);
	}
	Vector2 operator+(const Vector2& u) const
	{
		return Vector2(x + u.x, y + u.y);
	}
	Vector2 operator+(float k) const
	{
		return Vector2(x + k, y + k);
	}
	Vector2 operator-(float k) const
	{
		return Vector2(x - k, y - k);
	}
	float operator[](int i) const
	{
		if (i == 0)
			return x;
		return y;
	}
	float& operator[](int i)
	{
		if (i == 0)
			return x;
		return y;
	}
	friend std::ostream& operator<< (std::ostream& stream, const Vector2& u);
	inline Vector3 ToVector3(float yy) const
	{
		return Vector3(x, yy, y);
	}
	inline float Max() const
	{
		return KarstNSim::Max(x, y);
	}
	inline float Min() const
	{
		return KarstNSim::Min(x, y);
	}
};
inline std::ostream& operator<<(std::ostream& stream, const Vector2& u)
{
	stream << "(" << u.x << ", " << u.y << ");";
	return stream;
}
inline float Dot(const Vector2& u, const Vector2& v)
{
	return u.x * v.x + u.y * v.y;
}
inline float magnitude(const Vector2& u)
{
	return sqrt(u.x * u.x + u.y * u.y);
}
inline float squaredmagnitude(const Vector2& u)
{
	return u.x * u.x + u.y * u.y;
}
inline Vector2 Normalize(const Vector2& v)
{
	float kk = 1.0f / magnitude(v);
	return v * kk;
}
inline Vector2 operator-(const Vector2& v)
{
	return Vector2(-v.x, -v.y);
}
inline bool operator>(const Vector2& u, const Vector2& v)
{
	return (u.x > v.x) && (u.y > v.y);
}
inline bool operator<(const Vector2& u, const Vector2& v)
{
	return (u.x < v.x) && (u.y < v.y);
}
inline bool operator>=(const Vector2& u, const Vector2& v)
{
	return (u.x >= v.x) && (u.y >= v.y);
}
inline bool operator<=(const Vector2& u, const Vector2& v)
{
	return (u.x <= v.x) && (u.y <= v.y);
}
inline Vector2 Abs(const Vector2& u)
{
	return Vector2(u[0] > 0.0 ? u[0] : -u[0], u[1] > 0.0 ? u[1] : -u[1]);
}

/*!
\brief A 3x3 matrix structure for linear algebra operations.

This structure represents a 3x3 matrix and provides functionality for constructing matrices,
computing the determinant, and accessing matrix elements using the () operator.

\details
- The matrix is stored as a flat array of 9 elements in row-major order.
- Includes constructors for zero-initialized, vector-based, and element-based initialization.

\sa Vector3
*/
struct Matrix3
{
public:
	float r[9];

	Matrix3();
	Matrix3(const Vector3& a, const Vector3& b, const Vector3& c);
	Matrix3(float a00, float a01, float a02, float a10, float a11, float a12, float a20, float a21, float a22);
	float Determinant() const;
	float& operator() (int, int);
	float operator() (int, int) const;
};

inline Matrix3::Matrix3()
{
	for (int i = 0; i < 9; i++)
		r[i] = 0.0f;
}

inline Matrix3::Matrix3(const Vector3& a, const Vector3& b, const Vector3& c)
{
	r[0] = a[0];
	r[1] = a[1];
	r[2] = a[2];

	r[3] = b[0];
	r[4] = b[1];
	r[5] = b[2];

	r[6] = c[0];
	r[7] = c[1];
	r[8] = c[2];
}

inline Matrix3::Matrix3(float a00, float a01, float a02, float a10, float a11, float a12, float a20, float a21, float a22)
{
	r[0] = a00;
	r[1] = a01;
	r[2] = a02;

	r[3] = a10;
	r[4] = a11;
	r[5] = a12;

	r[6] = a20;
	r[7] = a21;
	r[8] = a22;
}

inline float Matrix3::Determinant() const
{
	return r[0] * r[4] * r[8] + r[1] * r[5] * r[6] + r[2] * r[3] * r[7] - r[2] * r[4] * r[6] - r[1] * r[3] * r[8] - r[0] * r[5] * r[7];
}

inline float& Matrix3::operator() (int i, int j)
{
	return r[i + j + j + j];
}

inline float Matrix3::operator() (int i, int j) const
{
	return r[i + j + j + j];
}


/*!
\brief A 4x4 matrix structure for linear algebra operations.

This structure represents a 4x4 matrix and provides functionality for constructing matrices,
computing the determinant, and accessing matrix elements using the () operator.

\details
- The matrix is stored as a flat array of 16 elements in row-major order.
- Includes a constructor for initializing from a 3x3 matrix.

\sa Matrix3
*/
struct Matrix4
{
public:
	float r[16];

	Matrix4();
	Matrix4(const Matrix3& m3);
	float& operator() (int, int);
	float operator() (int, int) const;
	float Determinant() const;
};

inline Matrix4::Matrix4()
{
	for (int i = 0; i < 16; i++)
		r[i] = 0.0f;
}

/*!
\brief Create an homogeneous matrix from a simple Matrix.
The translation and shear coefficients are set to 0.0.
\param a Matrix.
*/
inline Matrix4::Matrix4(const Matrix3& a)
{
	// Rotation and scale
	r[0] = a.r[0];
	r[1] = a.r[1];
	r[2] = a.r[2];
	r[4] = a.r[3];
	r[5] = a.r[4];
	r[6] = a.r[5];
	r[8] = a.r[6];
	r[9] = a.r[7];
	r[10] = a.r[8];

	// Translation
	r[3] = r[7] = r[11] = 0.0;

	// Shear
	r[12] = r[13] = r[14] = 0.0;

	// Scale
	r[15] = 1.0;
}

inline float& Matrix4::operator() (int i, int j)
{
	return r[i + j + j + j];
}

inline float Matrix4::operator() (int i, int j) const
{
	return r[i + j + j + j];
}

inline float Matrix4::Determinant() const
{
	const Matrix4& M = *this;
	return M(0, 0) * Matrix3(M(1, 1), M(1, 2), M(1, 3), M(2, 1), M(2, 2), M(2, 3), M(3, 1), M(3, 2), M(3, 3)).Determinant()
		- M(1, 0) * Matrix3(M(0, 1), M(0, 2), M(0, 3), M(2, 1), M(2, 2), M(2, 3), M(3, 1), M(3, 2), M(3, 3)).Determinant()
		+ M(2, 0) * Matrix3(M(0, 1), M(0, 2), M(0, 3), M(1, 1), M(1, 2), M(1, 3), M(3, 1), M(3, 2), M(3, 3)).Determinant()
		- M(3, 0) * Matrix3(M(0, 1), M(0, 2), M(0, 3), M(1, 1), M(1, 2), M(1, 3), M(2, 1), M(2, 2), M(2, 3)).Determinant();
}

/*!
\brief A structure representing a neighboring element with an index and distance.

This structure is used to store information about a neighbor, including its index,
a distance value, and provides comparison operators for sorting based on the distance.

\details
- The structure is primarily used for scenarios involving proximity-based queries,
  such as nearest-neighbor searches.
- The distance (`d`) can represent any metric, such as Euclidean distance.

\sa operator<, operator>
*/
struct Neighbour {
	int i;    /*!< Index of the neighbor. */
	float d;  /*!< Distance value representing proximity to a reference point. */
	bool operator<(const Neighbour& nei) const {
		return d < nei.d;
	}
	bool operator>(const Neighbour& nei) const {
		return d > nei.d;
	}
};

/*!
\brief Removes specified neighbors from a sample set.

This function removes elements from the `s` vector (representing sample points) based on the indices
provided in the `neighbors` vector. Only valid indices within the range of `s` are considered.

\param[in,out] s A vector of `Vector3` representing the sample points. This vector is modified in-place.
\param[in] neighbors A vector of `Neighbour` objects specifying the indices of elements to be removed.

\details
- The function uses an `unordered_set` to efficiently store indices of points to remove.
- Points that are not listed in the indices are preserved and copied to a new vector.
- The `s` vector is updated to reflect the filtered points.

\warning Invalid indices (e.g., negative or out-of-bounds) in the `neighbors` vector are ignored.

\sa Neighbour
*/
inline void remove_neighbors(std::vector<Vector3>& s, const std::vector<Neighbour>& neighbors) {

	// Use an unordered_set to store unique indices to be removed
	std::unordered_set<int> indicesToRemove;

	// Collect the indices from the neighbors vector
	for (const auto& neighbor : neighbors) {
		if (neighbor.i >= 0 && neighbor.i < s.size()) {
			indicesToRemove.insert(neighbor.i);  // Ensure indices are valid
		}
	}

	// Create a new vector to store the points that we want to keep
	std::vector<Vector3> newSamples;

	// Copy the points that are not in indicesToRemove
	for (int i = 0; i < s.size(); ++i) {
		if (indicesToRemove.find(i) == indicesToRemove.end()) {
			newSamples.push_back(s[i]);
		}
	}

	// Replace the old samples vector with the new one
	s = std::move(newSamples);
}

/*!
\brief Represents a point cloud with KD-tree capabilities.

The `PointCloud` struct manages a collection of 3D points (`Vector3`) and provides
spatial query functionality using a KD-tree for nearest neighbor searches.

\details
- The KD-tree is built using the `nanoflann` library, which optimizes spatial queries.
- The struct supports operations such as finding nearest neighbors in 3D space.

\sa nanoflann::KDTreeSingleIndexAdaptor
*/
struct PointCloud {

public:
	std::vector<Vector3> cloud; /*!< The collection of points in the point cloud. */
	/*!
	\brief Alias for the KD-tree type used to index the point cloud.
	\details The KD-tree uses the L2 (Euclidean) distance metric for spatial queries.
	*/
	using my_kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<float, PointCloud>, PointCloud, 3>;
	std::unique_ptr<my_kd_tree_t> index; /*!< Pointer to the KD-tree index for the point cloud. */
	static const std::vector<Vector3> defaultCloud; /*!< Default point cloud data, if required. */


	PointCloud();
	PointCloud(const std::vector<Vector3>& cloud);
	//PointCloud(const PointCloud &) = default;
	PointCloud(const PointCloud& other);
	size_t kdtree_get_point_count() const;
	float kdtree_get_pt(const size_t idx, const size_t dim) const;
	template <class BBOX>
	bool kdtree_get_bbox(BBOX&) const;
	/*!
	\brief Finds the nearest neighbors of a query point within a specified range.

	\param[in] queryPoint The point to search neighbors for.
	\param[in] querryidx Index of the query point.
	\param[in] N The maximum number of neighbors to find.
	\param[in] distmax The maximum allowable distance for a neighbor to be considered.
	\param[out] result A vector of `Neighbour` objects containing the results of the query.

	\details
	- The function uses the KD-tree index to efficiently perform spatial queries.
	- Neighbors are sorted by distance in ascending order.

	\sa Neighbour
	*/
	void findNearestNeighbors(const Vector3& queryPoint, int querryidx, int N, float distmax, std::vector<Neighbour>& result) const;
};

// Constructor for PointCloud
inline PointCloud::PointCloud() :cloud(defaultCloud)
{
}

// Constructor for PointCloud
inline PointCloud::PointCloud(const std::vector<Vector3>& cloud)
	: cloud(cloud),
	index(nullptr)
{
	// Create the index outside of the initialization list
	index = std::make_unique<my_kd_tree_t>(3, *this, nanoflann::KDTreeSingleIndexAdaptorParams(10));
	//index->buildIndex();
}

inline PointCloud::PointCloud(const PointCloud& other)
	: cloud(other.cloud),
	index(nullptr)
{
	// Create the index outside of the initialization list
	index = std::make_unique<my_kd_tree_t>(3, *this, nanoflann::KDTreeSingleIndexAdaptorParams(10));
	//index->buildIndex(); // Uncomment if building the index is necessary after copying
}

// Must return the number of data points
inline size_t PointCloud::kdtree_get_point_count() const { return cloud.size(); }

// Returns the dim'th component of the idx'th point in the class:
// Since this is inlined and the "dim" argument is typically an immediate
// value, the
//  "if/else's" are actually solved at compile time.
inline float PointCloud::kdtree_get_pt(const size_t idx, const size_t dim) const
{
	if (dim == 0)
		return cloud[idx].x;
	else if (dim == 1)
		return cloud[idx].y;
	else
		return cloud[idx].z;
}

// Optional bounding-box computation: return false to default to a standard
// bbox computation loop.
//   Return true if the BBOX was already computed by the class and returned
//   in "bb" so it can be avoided to redo it again. Look at bb.size() to
//   find out the expected dimensionality (e.g. 2 or 3 for point clouds)
template <class BBOX>
inline bool PointCloud::kdtree_get_bbox(BBOX& /* bb */) const
{
	return false;
}

// Function to find nearest neighbors
// queryPoint : point from which neighbors should be found
// querryidx : index of the querryPoint inside the PointCloud (used to avoid considering that the point is neighbor of itself)
// N : max number of neighbors found (if unsure : use a high value)
// distmax : max distance from the querryPoint to find neighbors
// result : resulting neighbors vector

inline void PointCloud::findNearestNeighbors(const Vector3& queryPoint, int querryidx, int N, float distmax, std::vector<Neighbour>& result) const {

	// Initialize the results structure
	std::vector<int> indices(N);
	std::vector<float> dists(N);

	// Search for the nearest neighbors
	nanoflann::KNNResultSet<float, int> resultSet(N);
	resultSet.init(&indices[0], &dists[0]);
	nanoflann::SearchParameters params;
	params.sorted = true; // Sort the results based on distance

	//index.buildIndex();
	index->findNeighbors(resultSet, &queryPoint.x, params);

	// Process the results
	for (size_t i = 0; i < indices.size(); ++i) {
		int index = indices[i];
		float distance = dists[i];

		// Check if the distance is within the specified range
		if (distance <= distmax && querryidx != index) {
			// Add the neighbor to the result vector
			result.push_back(Neighbour{ index, distance });
		}
	}
}

/*!
\brief A 2D array implementation with contiguous memory storage for high performance.

The `Array2D` class provides a 2D array structure with contiguous memory storage, making it more efficient for certain operations compared to a `std::vector<std::vector<T>>`.

\details
- The class is particularly suited for large simulations where the array may need frequent resizing or erasure.
- Its contiguous memory layout ensures faster operations like erasing or overwriting the entire data compared to nested vectors.
- Efficient resizing and memory management can significantly reduce computation time for large datasets.

\warning
- Erasing the `Array2D` is much faster than clearing a vector of vectors because the memory is stored contiguously, leading to reduced fragmentation and overhead.

\tparam T The type of elements stored in the array.
*/
template<typename T>
class Array2D {
public:
	// Default constructor
	Array2D() : x_dim(0), y_dim(0) {}

	// Constructor with dimensions and default value
	Array2D(size_t x, size_t y, const T& default_value = T())
		: x_dim(x), y_dim(y), data(x * y, default_value) {}

	// Constructor with multiple vectors
	Array2D(const std::vector<T>& vec1, const std::vector<T>& vec2) {
		// Check if all vectors have the same size
		if (vec1.size() != vec2.size()) {
			throw std::invalid_argument("All input vectors must have the same size.");
		}

		x_dim = 2; // Two rows
		y_dim = vec1.size();
		data.resize(x_dim * y_dim);

		// Copy vec1 and vec2 into the data vector
		std::copy(vec1.begin(), vec1.end(), data.begin());
		std::copy(vec2.begin(), vec2.end(), data.begin() + y_dim);
	}

	// Copy constructor
	Array2D(const Array2D& other)
		: x_dim(other.x_dim), y_dim(other.y_dim), data(other.data) {}

	// Copy assignment operator
	Array2D& operator=(const Array2D& other) {
		if (this != &other) {
			x_dim = other.x_dim;
			y_dim = other.y_dim;
			data = other.data;
		}
		return *this;
	}

	// Move constructor
	Array2D(Array2D&& other) noexcept
		: x_dim(other.x_dim), y_dim(other.y_dim), data(std::move(other.data)) {
		other.x_dim = 0;
		other.y_dim = 0;
	}

	// Move assignment operator
	Array2D& operator=(Array2D&& other) noexcept {
		if (this != &other) {
			x_dim = other.x_dim;
			y_dim = other.y_dim;
			data = std::move(other.data);
			other.x_dim = 0;
			other.y_dim = 0;
		}
		return *this;
	}

	// Resize method to change the number of rows while keeping the same number of columns
	void resize(size_t new_x_dim) {
		if (new_x_dim == x_dim) return; // No change needed

		// Resize the data vector
		std::vector<T> new_data(new_x_dim * y_dim, T());

		// Copy old data to new data vector
		for (size_t i = 0; i < std::min(new_x_dim, x_dim); ++i) {
			std::copy(data.begin() + i * y_dim, data.begin() + (i + 1) * y_dim, new_data.begin() + i * y_dim);
		}

		// Update dimensions and data
		x_dim = new_x_dim;
		data = std::move(new_data);
	}

	// Resize method
	void resize(size_t new_x_dim, size_t new_y_dim, const T& default_value = T()) {
		std::vector<T> new_data(new_x_dim * new_y_dim, default_value);

		size_t min_x_dim = std::min(new_x_dim, x_dim);
		size_t min_y_dim = std::min(new_y_dim, y_dim);

		for (size_t x = 0; x < min_x_dim; ++x) {
			for (size_t y = 0; y < min_y_dim; ++y) {
				new_data[x * new_y_dim + y] = data[x * y_dim + y];
			}
		}

		x_dim = new_x_dim;
		y_dim = new_y_dim;
		data = std::move(new_data);
	}

	// Get number of rows
	size_t size() const {
		return x_dim;
	}

	// Get number of columns
	size_t cols() const {
		return y_dim;
	}

	// Access entire row as a vector
	std::vector<T> row(size_t x) const {
		if (x >= x_dim) {
			throw std::out_of_range("Index out of range");
		}
		return std::vector<T>(data.begin() + x * y_dim, data.begin() + (x + 1) * y_dim);
	}

	// Access element
	T& at(size_t x, size_t y) {
		if (x >= x_dim || y >= y_dim) {
			throw std::out_of_range("Index out of range");
		}
		return data[x * y_dim + y];
	}

	const T& at(size_t x, size_t y) const {
		if (x >= x_dim || y >= y_dim) {
			throw std::out_of_range("Index out of range");
		}
		return data[x * y_dim + y];
	}

	// Access element using operator()
	T& operator()(size_t x, size_t y) {
		return data[x * y_dim + y];
	}

	const T& operator()(size_t x, size_t y) const {
		return data[x * y_dim + y];
	}

	/*!
	\brief A helper class for accessing and manipulating rows within the `Array2D`.

	The `Row` class provides a convenient interface for accessing, modifying, and iterating through individual rows of the 2D array.
	It operates directly on the underlying memory of the `Array2D` for efficiency.

	\details
	- Provides bounds-checked access methods (`at`) to ensure safe element access.
	- Supports iterator-based traversal and assignment of rows.
	- Enables swapping or clearing rows efficiently without additional memory allocation.
	- Ensures that all row operations respect the contiguous memory layout of the parent `Array2D`.

	\tparam T The type of elements stored in the row.
	*/
	class Row {
	public:
		explicit Row(T* start, size_t length) : start(start), length(length) {}

		// Access element with bounds checking
		T& at(size_t y) {
			if (y >= length) {
				throw std::out_of_range("Index out of range");
			}
			return start[y];
		}

		const T& at(size_t y) const {
			if (y >= length) {
				throw std::out_of_range("Index out of range");
			}
			return start[y];
		}

		T& operator[](size_t y) {
			if (y >= length) {
				throw std::out_of_range("Index out of range");
			}
			return start[y];
		}

		const T& operator[](size_t y) const {
			if (y >= length) {
				throw std::out_of_range("Index out of range");
			}
			return start[y];
		}

		size_t size() const {
			return length;
		}

		// Begin iterator for the row
		T* begin() {
			return start;
		}

		// End iterator for the row
		T* end() {
			return start + length;
		}

		// Const begin iterator for the row
		const T* begin() const {
			return start;
		}

		// Const end iterator for the row
		const T* end() const {
			return start + length;
		}

		// Swap function
		void swap(Row& other) {
			if (this->length != other.length) {
				throw std::logic_error("Cannot swap rows of different lengths");
			}
			std::swap_ranges(start, start + length, other.start);
		}

		void assign(const std::vector<T>& vec) {
			if (vec.size() != length) {
				throw std::invalid_argument("Vector size does not match row size.");
			}
			std::copy(vec.begin(), vec.end(), start);
		}

		// Check if the row is empty
		bool empty() const {
			return length == 0;
		}

		// Clear the row
		void clear() {
			std::fill(start, start + length, T());
		}

	public:
		T* start;
		size_t length;
	};

	// Access row as Row
	Row operator[](size_t x) {
		if (x >= x_dim) {
			throw std::out_of_range("Index out of range");
		}
		return Row(&data[x * y_dim], y_dim);
	}

	const Row operator[](size_t x) const {
		if (x >= x_dim) {
			throw std::out_of_range("Index out of range");
		}
		return Row(&data[x * y_dim], y_dim);
	}

	// Swap rows
	void swap_rows(size_t row1, size_t row2) {
		if (row1 >= x_dim || row2 >= x_dim) {
			throw std::out_of_range("Row index out of range");
		}
		if (row1 != row2) {
			(*this)[row1].swap((*this)[row2]);
		}
	}

	void set_row(size_t x, const std::vector<T>& vec) {
		if (x >= x_dim) {
			throw std::out_of_range("Row index out of range");
		}
		(*this)[x].assign(vec);
	}

	// Check if the entire 2D array is empty
	bool empty() const {
		return x_dim == 0 || y_dim == 0;
	}

public:
	size_t x_dim, y_dim;
	std::vector<T> data;
};

/*!
\brief Sets flags in a 2D array based on indices provided for each flag.

\details
This function processes a list of sample points and a corresponding list of flag indices to populate a 2D array of flags.
Each flag corresponds to a specific condition, and the indices for which the condition is true are stored in `on_wt_flags_idx`.
The function initializes all flags to `false` and sets them to `true` for the indices specified in `on_wt_flags_idx`.

\param[in] samples The list of sample points. (Not used directly but ensures consistency in dimensions.)
\param[in] on_wt_flags_idx A vector where each element is a list of indices that correspond to samples for which a specific flag should be set to `true`.
\param[out] samples_on_wt_flags A 2D array of flags (`true`/`false`), where rows correspond to samples and columns to different flags.

*/
inline void get_flags_from_indices(const std::vector<Vector3>& samples, const std::vector<std::vector<int>>& on_wt_flags_idx, Array2D<char>& samples_on_wt_flags) {

	samples_on_wt_flags.resize(samples.size(), on_wt_flags_idx.size(), false); // Initialize all flags to false

   // Set flags to true for indices present in on_wt_flags_idx
	for (size_t i = 0; i < on_wt_flags_idx.size(); ++i) {
		for (int idx : on_wt_flags_idx[i]) {
			if (idx >= 0 && idx < samples.size()) {
				samples_on_wt_flags(idx, i) = true;
			}
		}
	}
}
