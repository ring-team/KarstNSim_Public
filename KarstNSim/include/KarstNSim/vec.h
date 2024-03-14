/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr for some generic functions and modifications to other classes
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

Copyright (c) 2021 Axel Paris
Author : Axel Paris, for classes Vector2i, Vector2, Vector3, Vector4
If you use this code, pleace cite : Paris et al., 2021, Computer Graphic Forum.

***************************************************************/

#pragma once

/**
@file vec.h
@brief Utility functions for simple vector objects.
@authors Axel PARIS, Augustin GOUY
**/

#include <KarstNSim/nanoflann.hpp>
#include <cmath>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <algorithm>

/* Forward Declarations */
struct Vector2i;
struct Vector2;
struct Vector3;
struct Vector4;

// Maths utility
namespace KarstNSim
{
	/*!
	\brief Clamps a value between two limits
	\param x Value
	\param a, b Interval values.
	\return Real.
	*/
	inline double Clamp(double x, double a = 0.0f, double b=1.0f)
	{
		return x < a ? a : x > b ? b : x;
	}

	/*!
	\brief Create a linear step.
	\param x Value
	\param a, b Interval values.
	\return Real in unit inverval.
	*/
	inline double Step(const double& x, const double& a, const double& b)
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
	inline T Lerp(T a, T b, double t)
	{
		return (a * (1.0f - t)) + (b * t);
	}

	inline double Abs(double a)
	{
		return a < 0 ? -a : a;
	}

	inline double CubicSmoothCompact(double x, double r)
	{
		return (x > r) ? 0.0f : (1.0f - x / r) * (1.0f - x / r) * (1.0f - x / r);
	}
	template <class T> T _square(T a) { return a * a; };

	inline double CubicSigmoid(double x, double r, double t)
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
				double y = -x;
				return  -(y * (1.0f + y * ((3.0f * t - 2.0f * r) / (r * r) + y * (r - 2.0f * t) / (r * r * r))));
			}
			else
			{
				return -t;
			}
		}
	}

	inline double CubicSmooth(double x)
	{
		return x * x * (3.0 - 2.0 * x);
	}

	inline double CubicSmooth( double x, double r)
	{
		return (1.0 - x / r) * (1.0 - x / r) * (1.0 - x / r);
	}

	inline double CubicSmoothStep(double x, double a, double b)
	{
		if (x < a)
			return 0.0;
		else if (x > b)
			return 1.0;
		else
			return 1.0 - CubicSmooth((x - a) * (x - a), (b - a) * (b - a));
	}

	inline double QuinticSmooth(double t)
	{
		return pow(t, 3.0) * (t * (t * 6.0 - 15.0) + 10.0);
	}

	inline double Pow(double x, double e)
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

	inline double calculateCircumradius(double a, double b, double c) {
		double s = (a + b + c) / 2;
		double circumradius = (a * b * c) / (4 * sqrt(s * (s - a) * (s - b) * (s - c)));
		return circumradius;
	}

	inline bool isTriangleValid(double a, double b, double c) {
		// Choose a small multiplier for epsilon (e.g., 1e-5)
		const double epsilonMultiplier = 1e-5;

		// Calculate a relative epsilon based on the magnitudes of the side lengths
		double epsilon = epsilonMultiplier * std::max({ std::abs(a), std::abs(b), std::abs(c) });

		// Check for non-negative side lengths
		if (a < 0 || b < 0 || c < 0) {
			return false;
		}

		// Check the Triangle Inequality with the relative epsilon
		if (a + b > c + epsilon && a + c > b + epsilon && b + c > a + epsilon) {
			// Check for degenerate cases (zero area) with the relative epsilon
			double s = (a + b + c) / 2;
			double area = sqrt(s * (s - a) * (s - b) * (s - c));
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
	double x, y, z;

	Vector3() : x(0.0f), y(0.0f), z(0.0f) { }
	Vector3(double n) : x(n), y(n), z(n) {}
	Vector3(double x, double y, double z) : x(x), y(y), z(z) {}
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
	Vector3 operator*= (double f)
	{
		return Vector3(x * f, y * f, z * f);
	}
	Vector3 operator/= (double f)
	{
		return Vector3(x / f, y / f, z / f);
	}
	Vector3 operator*(const Vector3& u) const
	{
		return Vector3(x * u.x, y * u.y, z * u.z);
	}
	Vector3 operator*(double k) const
	{
		return Vector3(x * k, y * k, z * k);
	}
	Vector3 operator/(double k) const
	{
		return Vector3(x / k, y / k, z / k);
	}
	bool operator==(const Vector3& u) const
	{
		double eps = 1e-3;
		return (abs(x - u.x)<eps && abs(y - u.y) < eps && abs(z - u.z) < eps);
	}
	bool operator!=(const Vector3& u) const
	{
		double eps = 1e-3;
		return (abs(x - u.x) > eps || abs(y - u.y) > eps || abs(z - u.z) > eps);
	}
	Vector3 operator-(const Vector3& u) const
	{
		return Vector3(x - u.x, y - u.y, z - u.z);
	}
	Vector3 operator+(const Vector3& u) const
	{
		return Vector3(x + u.x, y + u.y, z + u.z);
	}
	Vector3 operator+(double k) const
	{
		return Vector3(x + k, y + k, z + k);
	}
	double operator[](int i) const
	{
		if (i == 0)
			return x;
		else if (i == 1)
			return y;
		return z;
	}
	double& operator[](int i)
	{
		if (i == 0)
			return x;
		else if (i == 1)
			return y;
		return z;
	}
	friend std::ostream& operator<<(std::ostream& stream, const Vector3& u);
	inline double Max() const
	{
		return KarstNSim::Max(KarstNSim::Max(x, y), z);
	}
	inline double Min() const
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

inline std::ostream& operator<<(std::ostream& stream, const Vector3& u)
{
	stream << "(" << u.x << ", " << u.y << ", " << u.z << ");";
	return stream;
}
inline Vector3 Cross(const Vector3& u, const Vector3& v)
{
	return Vector3((u.y * v.z) - (u.z * v.y), (u.z * v.x) - (u.x * v.z), (u.x * v.y) - (u.y * v.x));
}

inline double Dot2D(const Vector3& a, const Vector3& b) {
	return a.x * b.x + a.y * b.y;
}

inline double Dot(const Vector3& u, const Vector3& v)
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

inline std::vector<std::vector<bool>> get_flags_from_indices(const std::vector<Vector3>& samples, const std::vector<std::vector<int>>& on_wt_flags_idx) {

	std::vector<std::vector<bool>> samples_on_wt_flags(samples.size(), std::vector<bool>(on_wt_flags_idx.size(), false)); // Initialize all flags to false

	// Set flags to true for indices present in on_wt_flags_idx
	for (size_t i = 0; i < on_wt_flags_idx.size(); ++i) {
		for (int idx : on_wt_flags_idx[i]) {
			if (idx >= 0 && idx < samples.size()) {
				samples_on_wt_flags[idx][i] = true;
			}
		}
	}

	return samples_on_wt_flags;
}

inline double magnitude2D(const Vector3& u)
{
	return sqrt(u.x * u.x + u.y * u.y);
}

inline double squaredmagnitude2D(const Vector3& u)
{
	return u.x * u.x + u.y * u.y;
}

inline double magnitude(const Vector3& u)
{
	return sqrt(u.x * u.x + u.y * u.y + u.z * u.z);
}

inline double squaredmagnitude(const Vector3& u)
{
	return u.x * u.x + u.y * u.y + u.z * u.z;
}

inline Vector3 mid(const Vector3& u, const Vector3& v)
{
	return Vector3((u.x + v.x) / 2, (u.y + v.y) / 2, (u.z + v.z) / 2);
}

inline void squaredistance_to_segment(const Vector3& pt, const Vector3& p1, const Vector3& p2, double& squaredDistance, Vector3& closestPoint)
{
	Vector3 v = p2 - p1;
	Vector3 w = pt - p1;
	double eps = 1e-5;

	double c1 = Dot(w, v);
	if (c1 <= eps) {
		squaredDistance = squaredmagnitude(w);
		closestPoint = p1;
		return;
	}

	double c2 = Dot(v, v);
	if (c2 <= c1) {
		squaredDistance = squaredmagnitude(pt - p2);
		closestPoint = p2;
		return;
	}

	double b = c1 / c2;
	Vector3 pb = p1 + v * b;
	squaredDistance = squaredmagnitude(pt - pb);
	closestPoint = pb;
}

inline void squaredistance_to_segment2D(const Vector3& pt, const Vector3& p1, const Vector3& p2, double& squaredDistance, Vector3& closestPoint)
{
	Vector3 v = p2 - p1;
	Vector3 w = pt - p1;
	double eps = 1e-5;

	double c1 = Dot2D(w, v);
	if (c1 <= eps) {
		squaredDistance = squaredmagnitude2D(w);
		closestPoint = p1;
		return;
	}

	double c2 = Dot2D(v, v);
	if (c2 <= c1) {
		squaredDistance = squaredmagnitude2D(pt - p2);
		closestPoint = p2;
		return;
	}

	double b = c1 / c2;
	Vector3 pb = p1 + v * b;
	squaredDistance = squaredmagnitude2D(pt - pb);
	closestPoint = pb;
}


inline double magnitudeNormalized(const Vector3& u, double du, double dv, double dw)
{
	return sqrt((u.x * u.x) / (du*du) + (u.y * u.y) / (dv*dv) + (u.z * u.z) / (dw*dw));
}

inline double squaredmagnitudeNormalized(const Vector3& u,double one_over_delta_x_squared, double one_over_delta_y_squared,double one_over_delta_z_squared)
{
	return (u.x * u.x)*one_over_delta_x_squared + (u.y * u.y)*one_over_delta_y_squared + (u.z * u.z)*one_over_delta_z_squared;
}

inline Vector3 Normalize(const Vector3& v)
{
	double magSquared = v.x * v.x + v.y * v.y + v.z * v.z;

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
inline Vector3 operator*(double a, const Vector3& v)
{
	return v * a;
}
inline Vector3 Abs(const Vector3& u)
{
	return Vector3(u[0] > 0.0 ? u[0] : -u[0], u[1] > 0.0 ? u[1] : -u[1], u[2] > 0.0 ? u[2] : -u[2]);
}

// Function to standardize values to a specified range
inline void standardize_to_range(std::vector<double>& vec, double min_value, double max_value) {
	// Find the minimum and maximum values in the vector
	double min_val = *std::min_element(vec.begin(), vec.end());
	double max_val = *std::max_element(vec.begin(), vec.end());

	// Compute the range of the vector values
	double range = max_val - min_val;

	// Iterate through each element in the vector
	for (int i = 0; i < vec.size(); ++i) {
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

//this method alows to reformat the (0,1)-range random values of a point in a flattened grid of dimensions dim1*dim2*dim3
// This operation is often called raveling.

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

inline int64_t flatten_indices2(const int& dim1, const int& dim2, const int& dim3, const double& dimx, const double& dimy, const double& dimz, const Vector3& point) {

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

// This method allows to get the indices in a 3D grid from a flattened array index
// This operation is often called unraveling. Note : you don't really need the 3rd dimension (the smallest one)

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

	return Vector3(i,j,k);
}

// This allows to get the equivalent index of a grid in a new grid, in a raveled context (1D).
//You can define if the indices start at 0 or 1 in option (0 by default)

inline int64_t convert_index(int dim1_in, int dim2_in, int dim3_in, int dim1_out, int dim2_out, int dim3_out, int64_t idx_in, int index_init_norm = 0) {
	Vector3 coord = unravel3D(idx_in, dim1_in, dim2_in, index_init_norm);
	int64_t i_out = static_cast<int64_t>(std::floor(coord.x * dim1_out / dim1_in));
	int64_t j_out = static_cast<int64_t>(std::floor(coord.y * dim2_out / dim2_in));
	int64_t k_out = static_cast<int64_t>(std::floor(coord.z * dim3_out / dim3_in));
	int64_t idx_out = i_out + j_out * dim1_out + k_out * dim1_out * dim2_out;
	return idx_out;
}

// this method finds the scanning area around a flattened index position, according to the local density r

inline std::vector<int64_t> scan_area(double r, double r_min, const Vector3& point, int dim) {
	int64_t flattened_index = flatten_indices(dim, dim, dim, point);
	Vector3 coord = unravel3D(flattened_index, dim, dim);
	std::vector<int64_t> scanned_indices;
	int region_size = static_cast<int>(std::ceil(2 * r * std::sqrt(3) / r_min)); // we're looking to scan all the elements at a distance less than region_size of the index
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

inline std::vector<int64_t> scan_area2(const double& r, const double& r_min, const Vector3& point, const int& dim1, const int& dim2,
	const int& dim3, const double& dimx, const double& dimy, const double& dimz) {

	int64_t flattened_index = flatten_indices2(dim1, dim2, dim3, dimx, dimy, dimz, point);
	Vector3 coord = unravel3D(flattened_index, dim1, dim2);
	std::vector<int64_t> scanned_indices;
	int region_size = static_cast<int>(std::ceil(2 * r * std::sqrt(3) / r_min)); // we're looking to scan all the elements at a distance less than region_size of the index
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

// Checks if point is inside triangle IN MAP VIEW (2D, so no z check!!)

inline bool isInsideTriangle(const Vector3& A, const Vector3& B, const Vector3& C, const Vector3& P, const double tolerance = 1e-2) {
	// Calculate barycentric coordinates
	const double detT = (B.y - C.y) * (A.x - C.x) + (C.x - B.x) * (A.y - C.y);
	const double alpha = ((B.y - C.y) * (P.x - C.x) + (C.x - B.x) * (P.y - C.y)) / detT;
	const double beta = ((C.y - A.y) * (P.x - C.x) + (A.x - C.x) * (P.y - C.y)) / detT;
	const double gamma = 1.0 - alpha - beta;

	// Check if the point is inside the triangle with tolerance
	return alpha >= -tolerance && beta >= -tolerance && gamma >= -tolerance;
}

// Function to push back a Vector3 if not already in the vector with a tolerance
inline void pushBackIfNotPresent(std::vector<Vector3>& vec, const Vector3& newVector, double tolerance = 1e-5) {
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

inline double scalar_product(const Vector3& v1, const Vector3& v2) {
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

// Check if a point is below a triangle
inline bool is_point_under_triangle(const Vector3& new_horizon_pt, const Vector3& P0, const Vector3& P1, const Vector3& P2, const double tolerance = 1) {
	Vector3 edge1 = { P1.x - P0.x, P1.y - P0.y, P1.z - P0.z };
	Vector3 edge2 = { P2.x - P0.x, P2.y - P0.y, P2.z - P0.z };
	Vector3 normal;
	vect_product(edge1, edge2,normal);
	double d = -scalar_product(normal, P0);
	double res = scalar_product(new_horizon_pt, normal) + d;
	if (std::abs(res) <= tolerance) {
		return true;  // Point is coplanar or very close
	}
	else if (res > tolerance) {
		return true;  // Point is below the triangle
	}
	else {
		return false;  // Point is above the triangle
	}
}


inline double barycentric(const Vector3& pt, const Vector3& A, const Vector3& B, const Vector3& C) {

	// Compute the normal vector of the plane defined by the triangle
	Vector3 normal = Cross(B - A, C - A);

	// Normalize the normal vector
	const double normalLength = std::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
	normal.x /= normalLength;
	normal.y /= normalLength;
	normal.z /= normalLength;

	// Compute the signed distance to the plane
	const double signedDistance = Dot(normal, pt - A);

	// Return the absolute value of the signed distance (orthogonal projection distance)
	return std::abs(signedDistance);
}


// Function to compute the signed distance from a point to the plane defined by a triangle
inline double signed_distance_to_plane(const Vector3& pt, const Vector3& A, const Vector3& B, const Vector3& C) {
	const Vector3 normal = Normalize(Cross((B - A),(C - A)));
	return Dot(normal,pt - A);
}

/* Vector2 */
struct Vector2
{
public:
	double x, y;

	explicit Vector2() : x(0.0f), y(0.0f) { }
	explicit Vector2(const Vector2i& v) : x(double(v.x)), y(double(v.y)) { }
	explicit Vector2(double n) : x(n), y(n) { }
	explicit Vector2(double x, double y) : x(x), y(y) { }
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
	Vector2 operator*= (double f)
	{
		return Vector2(x * f, y * f);
	}
	Vector2 operator/= (double f)
	{
		return Vector2(x / f, y / f);
	}
	Vector2 operator*(const Vector2& v) const
	{
		return Vector2(x * v.x, y * v.y);
	}
	Vector2 operator*(double k) const
	{
		return Vector2(x * k, y * k);
	}
	Vector2 operator/(double k) const
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
	Vector2 operator+(double k) const
	{
		return Vector2(x + k, y + k);
	}
	Vector2 operator-(double k) const
	{
		return Vector2(x - k, y - k);
	}
	double operator[](int i) const
	{
		if (i == 0)
			return x;
		return y;
	}
	double& operator[](int i)
	{
		if (i == 0)
			return x;
		return y;
	}
	friend std::ostream& operator<< (std::ostream& stream, const Vector2& u);
	inline Vector3 ToVector3(double yy) const
	{
		return Vector3(x, yy, y);
	}
	inline double Max() const
	{
		return KarstNSim::Max(x, y);
	}
	inline double Min() const
	{
		return KarstNSim::Min(x, y);
	}
};
inline std::ostream& operator<<(std::ostream& stream, const Vector2& u)
{
	stream << "(" << u.x << ", " << u.y << ");";
	return stream;
}
inline double Dot(const Vector2& u, const Vector2& v)
{
	return u.x * v.x + u.y * v.y;
}
inline double magnitude(const Vector2& u)
{
	return sqrt(u.x * u.x + u.y * u.y);
}
inline double squaredmagnitude(const Vector2& u)
{
	return u.x * u.x + u.y * u.y;
}
inline Vector2 Normalize(const Vector2& v)
{
	double kk = 1.0f / magnitude(v);
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

/* Matrix3 */
struct Matrix3
{
public:
	double r[9];

	Matrix3();
	Matrix3(const Vector3& a, const Vector3& b, const Vector3& c);
	Matrix3(double a00, double a01, double a02, double a10, double a11, double a12, double a20, double a21, double a22);
	double Determinant() const;
	double& operator() (int, int);
	double operator() (int, int) const;
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

inline Matrix3::Matrix3(double a00, double a01, double a02, double a10, double a11, double a12, double a20, double a21, double a22)
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

inline double Matrix3::Determinant() const
{
	return r[0] * r[4] * r[8] + r[1] * r[5] * r[6] + r[2] * r[3] * r[7] - r[2] * r[4] * r[6] - r[1] * r[3] * r[8] - r[0] * r[5] * r[7];
}

inline double& Matrix3::operator() (int i, int j)
{
	return r[i + j + j + j];
}

inline double Matrix3::operator() (int i, int j) const
{
	return r[i + j + j + j];
}


/* Matrix4 */
struct Matrix4
{
public:
	double r[16];

	Matrix4();
	Matrix4(const Matrix3& m3);
	double& operator() (int, int);
	double operator() (int, int) const;
	double Determinant() const;
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

inline double& Matrix4::operator() (int i, int j)
{
	return r[i + j + j + j];
}

inline double Matrix4::operator() (int i, int j) const
{
	return r[i + j + j + j];
}

inline double Matrix4::Determinant() const
{
	const Matrix4& M = *this;
	return M(0, 0) * Matrix3(M(1, 1), M(1, 2), M(1, 3), M(2, 1), M(2, 2), M(2, 3), M(3, 1), M(3, 2), M(3, 3)).Determinant()
		- M(1, 0) * Matrix3(M(0, 1), M(0, 2), M(0, 3), M(2, 1), M(2, 2), M(2, 3), M(3, 1), M(3, 2), M(3, 3)).Determinant()
		+ M(2, 0) * Matrix3(M(0, 1), M(0, 2), M(0, 3), M(1, 1), M(1, 2), M(1, 3), M(3, 1), M(3, 2), M(3, 3)).Determinant()
		- M(3, 0) * Matrix3(M(0, 1), M(0, 2), M(0, 3), M(1, 1), M(1, 2), M(1, 3), M(2, 1), M(2, 2), M(2, 3)).Determinant();
}

struct Neighbour {
	int i;
	double d;
	bool operator<(const Neighbour& nei) const {
		return d < nei.d;
	}
	bool operator>(const Neighbour& nei) const {
		return d > nei.d;
	}
};

struct PointCloud {

public:
		std::vector<Vector3> cloud;
		using my_kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud>, PointCloud, 3>;
		std::unique_ptr<my_kd_tree_t> index;
		static const std::vector<Vector3> defaultCloud;

		PointCloud();
		PointCloud(const std::vector<Vector3>& cloud);
		//PointCloud(const PointCloud &) = default;
		PointCloud(const PointCloud& other);
		size_t kdtree_get_point_count() const;
		double kdtree_get_pt(const size_t idx, const size_t dim) const;
		template <class BBOX>
		bool kdtree_get_bbox(BBOX&) const;
		void findNearestNeighbors(const Vector3& queryPoint, int querryidx, int N, double distmax, std::vector<Neighbour>& result) const;
};

// Constructor for PointCloud
inline PointCloud::PointCloud():cloud(defaultCloud)
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
inline double PointCloud::kdtree_get_pt(const size_t idx, const size_t dim) const
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
inline void PointCloud::findNearestNeighbors(const Vector3& queryPoint, int querryidx, int N, double distmax, std::vector<Neighbour>& result) const {
	// Initialize the results structure
	std::vector<int> indices(N);
	std::vector<double> dists(N);

	// Search for the nearest neighbors
	nanoflann::KNNResultSet<double, int> resultSet(N);
	resultSet.init(&indices[0], &dists[0]);
	nanoflann::SearchParameters params;
	params.sorted = true; // Sort the results based on distance

	//index.buildIndex();
	index->findNeighbors(resultSet, &queryPoint.x, params);

	// Process the results
	for (size_t i = 0; i < indices.size(); ++i) {
		int index = indices[i];
		double distance = dists[i];

		// Check if the distance is within the specified range
		if (distance <= distmax && querryidx != index) {
			// Add the neighbor to the result vector
			result.push_back(Neighbour{ index, distance });
		}
	}
}