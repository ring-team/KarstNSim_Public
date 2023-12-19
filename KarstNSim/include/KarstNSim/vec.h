/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr for some generic functions and modifications to other classes
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

Copyright (c) 2021 Axel Paris
Author : Axel Paris, for classes Vector2i, Vector2, Vector3, Vector4
If you use this code, pleace cite : Paris et al., 2021, Computer Graphic Forum.
Find the corresponding code at https://github.com/aparis69/Karst-Synthesis

***************************************************************/

#pragma once

/**
@file vec.h
@brief Utility functions for simple vector objects.
@authors Axel PARIS, Augustin GOUY
**/

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
	inline double Clamp(double x, double a = 0.0f, double b = 1.0f)
	{
		return x < a ? a : x > b ? b : x;
	}

	/*!
	\brief Create a linear step.
	\param x Value
	\param a, b Interval values.
	\return Real in unit inverval.
	*/
	inline double Step(double x, double a, double b)
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
		return x * x * (3.0f - 2.0f * x);
	}

	inline double CubicSmooth(double x, double r)
	{
		return (1.0f - x / r) * (1.0f - x / r) * (1.0f - x / r);
	}

	inline double CubicSmoothStep(double x, double a, double b)
	{
		if (x < a)
			return 0.0f;
		else if (x > b)
			return 1.0f;
		else
			return 1.0f - CubicSmooth((x - a) * (x - a), (b - a) * (b - a));
	}

	inline double QuinticSmooth(double t)
	{
		return pow(t, 3.0f) * (t * (t * 6.0f - 15.0f) + 10.0f);
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
inline double Dot(const Vector3& u, const Vector3& v)
{
	return u.x * v.x + u.y * v.y + u.z * v.z;
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
	double kk = 1.0f / magnitude(v);
	return v * kk;
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


//this method alows to reformat the (0,1)-range random values of a point in a flattened grid of dimensions dim1*dim2*dim3
// This operation is often called raveling.

inline int flatten_indices(int dim1, int dim2, int dim3, const Vector3& point) {
	int flattened_index = 0;
	int x = ((int)(ceil(point.x*dim1))) - 1; // values are ceiled (rounded up)
	int y = ((int)(ceil(point.y*dim2))) - 1;
	int z = ((int)(ceil(point.z*dim3))) - 1;
	flattened_index = x + y * dim1 + z * dim1 * dim2;
	return flattened_index;
}

inline int flatten_indices2(int dim1, int dim2, int dim3, double dimx, double dimy, double dimz, const Vector3& point) {
	int flattened_index = 0;
	int x = ((int)(ceil(point.x*dim1 / dimx))) - 1; // values are ceiled (rounded up)
	int y = ((int)(ceil(point.y*dim2 / dimy))) - 1;
	int z = ((int)(ceil(point.z*dim3 / dimz))) - 1;
	flattened_index = x + y * dim1 + z * dim1 * dim2;
	return flattened_index;
}

// This method allows to get the indices in a 3D grid from a flattened array index
// This operation is often called unraveling. Note : you don't really need the 3rd dimension (the smallest one)

inline Vector3 unravel3D(int idx_in, int intdim1, int intdim2, int index_init_norm = 0) {
	int i, j, k;
	if (index_init_norm == 0) {
		k = (int)(idx_in / (intdim1*intdim2));
		j = (int)((idx_in - k * intdim1*intdim2) / intdim1);
		i = (int)(idx_in - k * intdim1*intdim2 - j * intdim1);
	}
	else {
		k = (int)(ceil(idx_in / (intdim1*intdim2)));
		j = (int)(ceil((idx_in - k * intdim1*intdim2) / intdim1));
		i = (int)(ceil(idx_in - k * intdim1*intdim2 - j * intdim1));
	}
	return 	Vector3(i, j, k);
}

// This allows to get the equivalent index of a grid in a new grid, in a raveled context (1D).
//You can define if the indices start at 0 or 1 in option (0 by default)

inline int convert_index(int dim1_in, int dim2_in, int dim3_in, int dim1_out, int dim2_out, int dim3_out, int idx_in, int index_init_norm = 0) {
	Vector3 coord = unravel3D(idx_in, dim1_in, dim2_in, index_init_norm);
	int i_out = int(coord.x*dim1_out / dim1_in);
	int j_out = int(coord.y*dim2_out / dim2_in);
	int k_out = int(coord.z*dim3_out / dim3_in);
	int idx_out = i_out + j_out * dim1_out + k_out * dim1_out * dim2_out;
	return idx_out;
}

// this method finds the scanning area around a flattened index position, according to the local density r

inline std::vector<int> scan_area(double r, double r_min, const Vector3& point, int dim) {
	int flattened_index = flatten_indices(dim, dim, dim, point);
	Vector3 coord = unravel3D(flattened_index, dim, dim);
	std::vector<int> scanned_indices;
	int region_size = (int)(ceil(2 * r*sqrt(3) / r_min)); // we're looking to scan all the elements at a distance less than region_size of the index
	int x_i, y_i, z_i;
	for (x_i = -region_size; (x_i <= region_size); x_i++) {
		if ((0 <= coord.x + x_i) && (coord.x + x_i < dim)) { // we also check that those elements aren't outside the grid
			for (y_i = -region_size; (y_i <= region_size); y_i++) {
				if ((0 <= coord.y + y_i) && (coord.y + y_i < dim)) {
					for (z_i = -region_size; (z_i <= region_size); z_i++) {
						if ((0 <= coord.z + z_i) && (coord.z + z_i < dim)) {
							scanned_indices.push_back(flattened_index + x_i + y_i * (dim)+z_i * (dim)*(dim));
						}
					}
				}
			}
		}
	}
	return scanned_indices;
}

// version with varying dimensions

inline std::vector<int> scan_area2(double r, double r_min, const Vector3& point, int dim1, int dim2, int dim3, double dimx, double dimy, double dimz) {
	int flattened_index = flatten_indices2(dim1, dim2, dim3, dimx, dimy, dimz, point);
	Vector3 coord = unravel3D(flattened_index, dim1, dim2);
	std::vector<int> scanned_indices;
	int region_size = (int)(ceil(2 * r*sqrt(3) / r_min)); // we're looking to scan all the elements at a distance less than region_size of the index
	int x_i, y_i, z_i;
	for (x_i = -region_size; (x_i <= region_size); x_i++) {
		if ((0 <= coord.x + x_i) && (coord.x + x_i < dim1)) { // we also check that those elements aren't outside the grid
			for (y_i = -region_size; (y_i <= region_size); y_i++) {
				if ((0 <= coord.y + y_i) && (coord.y + y_i < dim2)) {
					for (z_i = -region_size; (z_i <= region_size); z_i++) {
						if ((0 <= coord.z + z_i) && (coord.z + z_i < dim3)) {
							scanned_indices.push_back(flattened_index + x_i + y_i * (dim1)+z_i * (dim1)*(dim2));
						}
					}
				}
			}
		}
	}
	return scanned_indices;
}

inline bool is_point_in_triangle(Vector3& new_horizon_pt, Vector3& A, Vector3& B, Vector3& origin) {
	double d1 = (B.x - A.x)*(new_horizon_pt.y - A.y) - (B.y - A.y)*(new_horizon_pt.x - A.x);
	double d2 = (B.x - A.x)*(origin.y - A.y) - (B.y - A.y)*(origin.x - A.x);
	if (d1 < 0 && d2 < 0) {
		return true;
	}
	else if (d1 > 0 && d2 > 0) {
		return true;
	}
	else {
		return false;
	}
}

inline void vect_product(Vector3 vec1, Vector3 vec2, Vector3& result) {
	result.x = vec1.y*vec2.z - vec1.z*vec2.y;
	result.y = vec1.z*vec2.x - vec1.x*vec2.z;
	result.z = vec1.x*vec2.y - vec1.y*vec2.x;
}

inline double scalar_product(Vector3 v1, Vector3 v2) {
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

inline bool is_point_under_triangle(Vector3& new_horizon_pt, Vector3& P0, Vector3& P1, Vector3& P2) {

	// preliminary test : if all 3 z values are above or below the point, no need to go further :
	if (P0.z >= new_horizon_pt.z && P1.z >= new_horizon_pt.z && P2.z >= new_horizon_pt.z) { // pz = z(triangle) is also accepted as under the triangle
		return true;
	}
	if (P0.z < new_horizon_pt.z && P1.z < new_horizon_pt.z && P2.z < new_horizon_pt.z) {
		return false;
	}

	Vector3 vec1(P1.x - P0.x, P1.y - P0.y, P1.z - P0.z);
	Vector3 vec2(P2.x - P0.x, P2.y - P0.y, P2.z - P0.z);
	Vector3 normal;
	vect_product(vec1, vec2, normal);
	Vector3 P0_new_horizon_pt(new_horizon_pt.x - P0.x, new_horizon_pt.y - P0.y, new_horizon_pt.z - P0.z);
	double res = scalar_product(P0_new_horizon_pt, normal);
	if (normal.z < 0) {//check orientation of the normal
		if (res >= 0) {
			return true;
		}
		else {
			return false;
		}
	}
	else if (res <= 0) {
		return true;
	}
	else {
		return false;
	}

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
