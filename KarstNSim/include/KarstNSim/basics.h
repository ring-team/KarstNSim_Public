/***************************************************************

Université de Lorraine - ANDRA - BRGM
Copyright(c) 2023 Université de Lorraine - ANDRA - BRGM. All Rights Reserved.
This code is published under the MIT License.
Author : Augustin Gouy - augustin.gouy@univ-lorraine.fr for classes Triangle, Surface, Box, Segment, Line
If you use this code, please cite : Gouy et al., 2024, Journal of Hydrology.

Copyright (c) 2021 Axel Paris 
Author : Axel Paris, for classes Random, Ray, Plane, Box2D, Sphere, ScalarField2D, Circle2D, Circle
If you use this code, pleace cite : Paris et al., 2021, Computer Graphic Forum.
Find the corresponding code at https://github.com/aparis69/Karst-Synthesis

***************************************************************/

#pragma once

/*!
\file basics.h
\brief Declaration of various classes representing objects used for (hydro)geological model and karst network representation.
\author Augustin GOUY, Axel PARIS
*/

#include "KarstNSim/vec.h"
#include <time.h>
#include <vector>
#include <algorithm>
#include <windows.h>

namespace KarstNSim {
// Random (Dirty, C-style)

class Random
{
public:
	/*!
	\brief Constructor.
	*/
	Random()
	{
		// Empty
	}

	/*!
	\brief Compute a random number in a given range.
	\param a min
	\param b max
	*/
	static inline double Uniform(double a, double b)
	{
		return a + (b - a) * Uniform();
	}

	/*!
	\brief Compute a uniform random number in [0, 1]
	*/
	static inline double Uniform()
	{
		return double(rand()) / RAND_MAX;
	}

	/*!
	\brief Compute a random positive integer.
	*/
	static inline int Integer()
	{
		return rand();
	}
};


// 3D Ray
class Ray
{
public:
	Vector3 o;
	Vector3 d;

	Ray();
	Ray(const Vector3& oo, const Vector3& dd);
	Vector3 operator()(double t) const;
};

inline Ray::Ray()
{

}

inline Ray::Ray(const Vector3& oo, const Vector3& dd)
{
	o = oo;
	d = dd;
}

inline Vector3 Ray::operator()(double t) const
{
	return o + t * d;
}


// Plane
class Plane
{
protected:
	Vector3 p;
	Vector3 n;
	double c;

public:
	Plane();
	Plane(const Vector3& p, const Vector3& n);

	Vector3 Normal() const;
	Vector3 Point() const;
	int Side(const Vector3& p) const;
	double Signed(const Vector3& p) const;
	static bool Intersection(const Plane& a, const Plane& b, const Plane& c, Vector3& p);
	static std::vector<Vector3> ConvexPoints(const std::vector<Plane>& planes);
};

/*!
\brief
*/
inline Plane::Plane()
{
	p = Vector3(0);
	n = Vector3(0);
	c = 0.0;
}

/*!
\brief Constructor from a point and a normal.
\param pp point on the plane
\param plane normal (should be normalized)
*/
inline Plane::Plane(const Vector3& pp, const Vector3& nn)
{
	p = pp;
	n = nn;
	c = Dot(p, n);
}

/*!
\brief Returns the normal vector of the plane.
*/
inline Vector3 Plane::Normal() const
{
	return n;
}

/*!
\brief Computes the signed distance to the plane.
\param p point
*/
inline double Plane::Signed(const Vector3& pp) const
{
	return Dot(n, pp) - c;
}

/*!
\brief Returns a point on the plane.
*/
inline Vector3 Plane::Point() const
{
	return p;
}

/*!
\brief Check if a point lies inside, outside or even on the plane.
The epsilon tolerance used for testing if the point is on the plane
is 10<SUP>-6</SUP>.
\param p Point.
*/
inline int Plane::Side(const Vector3& pp) const
{
	double r = Dot(n, pp) - c;

	// Epsilon test
	if (r > 1e-6)
		return 1;
	else if (r < -1e-6)
		return -1;
	else
		return 0;
}

/*!
\brief
*/
inline bool Plane::Intersection(const Plane& a, const Plane& b, const Plane& c, Vector3& p)
{
	double e = Matrix4(Matrix3(a.Normal(), b.Normal(), c.Normal())).Determinant();
	if (e < 1e-06)
		return false;
	p = (Dot(a.Point(), a.Normal()) * (Cross(b.Normal(), c.Normal()))) +
		(Dot(b.Point(), b.Normal()) * (Cross(c.Normal(), a.Normal()))) +
		(Dot(c.Point(), c.Normal()) * (Cross(a.Normal(), b.Normal())));
	p = p / (-e);
	return true;
}

/*!
\brief Compute the intersection of a set of half spaces.
This is a O(n^4) algorithm.
A better version in O(n^3) could be implemented: see Finding the Intersection
of Half-Spaces in Time O(n ln n). Preparata and Muller, Theoretical Computer Science 8, 45-55, 1979.
Returns a set of point representing the minimal convex polygon embedding all half spaces.
The function doesn't check if the intersection is bounded or not.
\param planes Set of planes.
*/
inline std::vector<Vector3> Plane::ConvexPoints(const std::vector<Plane>& planes)
{
	std::vector<Vector3> pts;
	for (int i = 0; i < int(planes.size()); i++)
	{
		for (int j = i + 1; j < int(planes.size()); j++)
		{
			for (int k = j + 1; k < int(planes.size()); k++)
			{
				Vector3 p;
				bool intersect = Intersection(planes[i], planes[j], planes[k], p);
				if (intersect)
				{
					bool isInside = true;
					for (int l = 0; l < int(planes.size()); l++)
					{
						// Do not check point ijk with one of its generating plane
						if (l == i || l == j || l == k)
							continue;
						int s = planes[l].Side(p);
						if (s > 0)
						{
							isInside = false;
							break;
						}
					}
					if (isInside)
						pts.push_back(p);
				}
			}
		}
	}
	return pts;
}


// Sphere. Spherical geometric element.
class Sphere
{
protected:
	Vector3 center;
	double radius;

public:
	Sphere(const Vector3& c, double r);

	double Distance(const Vector3& p) const;
	bool Contains(const Vector3& p) const;
	Vector3 RandomInside() const;
	Vector3 RandomSurface() const;
	Vector3 Center() const;
	double Radius() const;
	void Poisson(std::vector<Vector3>& p, double r, int n) const;
	void Poisson_adapted(std::vector<Vector3>& p, double r, int n, std::vector<Sphere>& spheres) const;
};

/*!
\brief Constructor.
\param c center
\param r radius
*/
inline Sphere::Sphere(const Vector3& c, double r) : center(c), radius(r)
{

}

/*!
\brief Compute the distance between a point and a sphere.
If the point lies inside the sphere, the distance is considered to be 0.
\param p world point
*/
inline double Sphere::Distance(const Vector3& p) const
{
	double a = Dot((p - center), (p - center));
	if (a < radius * radius)
	{
		return 0.0;
	}
	a = sqrt(a) - radius;
	a *= a;
	return a;
}

/*!
\brief Returns true of the point lies inside the sphere (ie. Distance(p) == 0)
\param p world point
*/
inline bool Sphere::Contains(const Vector3& p) const
{
	return Distance(p) == 0;
}

// 3D Triangle
class Triangle
{
private:
	int pts[3];

public:
	Triangle();
	Triangle(const int& a, const int& b, const int& c);
	//Vector3 Center() const;
	//Vector3 Normal() const;
	int point(int i) const;
	//Vector3 Closest(Vector3 pt) const;
};

/*!
\brief Default Constructor.
*/
inline Triangle::Triangle()
{
}

/*!
\brief Constructor from three points.
\param a, b, c points.
*/
inline Triangle::Triangle(const int& a, const int& b, const int& c)
{
	pts[0] = a;
	pts[1] = b;
	pts[2] = c;
}

///*!
//\brief Computes and returns the center of the triangle.
//*/
//inline Vector3 Triangle::Center() const
//{
//	return (pts[0] + pts[1] + pts[2]) / 3.0f;
//}
//
///*!
//\brief Computes and returns the normal of the triangle.
//*/
//inline Vector3 Triangle::Normal() const
//{
//	return Normalize(Cross(pts[1] - pts[0], pts[2] - pts[0]));
//}

/*!
\brief Returns one of the point of the triangle.
\param i point index
*/
inline int Triangle::point(int i) const
{
	return pts[i];
}

//inline Vector3 Triangle::Closest(Vector3 pt) const
//{
//	double dist1 = (pt.x-pts[0].x)*(pt.x - pts[0].x)+(pt.y - pts[0].y)*(pt.y - pts[0].y)+(pt.z - pts[0].z)*(pt.z - pts[0].z);
//	double dist2 = (pt.x - pts[1].x)*(pt.x - pts[1].x) + (pt.y - pts[1].y)*(pt.y - pts[1].y) + (pt.z - pts[1].z)*(pt.z - pts[1].z);
//	double dist3 = (pt.x - pts[2].x)*(pt.x - pts[2].x) + (pt.y - pts[2].y)*(pt.y - pts[2].y) + (pt.z - pts[2].z)*(pt.z - pts[2].z);
//	if (dist1<dist2 && dist1<dist3) return pt[0]
//	else if (dist2<dist3) return pt[1]
//	else return pt[2]
//}

class Surface
{
private:
	std::vector<Vector3> pts_;
	std::vector<Triangle> trgls_;
	Vector3 b_min_;
	Vector3 b_max_;
	double largest_edge_length_;

public:
	Surface();
	Surface(std::vector<Vector3>, std::vector<Triangle>);
	int get_nb_trgls();
	int get_nb_pts();
	Triangle get_triangle(int i);
	const Vector3& get_node(int i) const {
		return pts_[i];
	};
	//Vector3 get_node(int i);
	Vector3 nearest(Vector3 p1, Vector3 p2, Vector3 ptest);
	Vector3 nearest(Vector3 p1, Vector3 p2, Vector3 p3, Vector3 ptest);
	Vector3 get_boundbox_min();
	Vector3 get_boundbox_max();
	double get_largest_edge_length();
	//Vector3 Center() const;
	//Vector3 Normal() const;
	//Vector3 Point(int i) const;

};

/*!
\brief Default Constructor.
*/
inline Surface::Surface()
{
}

/*!
\brief Constructor from a vector of triangles.
\param trgls, vector of triangles
*/
inline Surface::Surface(std::vector<Vector3> pts, std::vector<Triangle> trgls)
{
	pts_ = pts;
	trgls_ = trgls;
	double x_min = pts[0].x;
	double x_max = pts[0].x;
	double y_min = pts[0].y;
	double y_max = pts[0].y;
	double z_min = pts[0].z;
	double z_max = pts[0].z;
	for (int i = 0 ; i < int(pts.size()); i++) {
		if (pts[i].x < x_min) x_min = pts[i].x;
		if (pts[i].x > x_max) x_max = pts[i].x;
		if (pts[i].y < y_min) y_min = pts[i].y;
		if (pts[i].y > y_max) y_max = pts[i].y;
		if (pts[i].z < z_min) z_min = pts[i].z;
		if (pts[i].z > z_max) z_max = pts[i].z;
	}
	b_min_ = Vector3(x_min, y_min, z_min); // this is the coordinates of the bottom corner in all 3 directions
	b_max_ = Vector3(x_max, y_max, z_max); // this is the coordinates of the top corner in all 3 directions

	// another parameter that is practical for TINs is the length of the largest circumradius among all its triangles (in map view)
	// We compute it directly when creating the surface.

	double largest_found = 0.;
	for (int trgl = 0; trgl < trgls.size(); trgl++) {
		Triangle trgl_i = trgls[trgl];
		int pt1 = trgl_i.point(0), pt2 = trgl_i.point(1), pt3 = trgl_i.point(2);
		Vector3 p1 = pts[pt1], p2 = pts[pt2], p3 = pts[pt3];
		// 2D dist between p1 and p2
		double dist1 = magnitude2D(p1-p2);
		// 2D dist between p2 and p3
		double dist2 = magnitude2D(p2 - p3);
		// 2D dist between p1 and p3
		double dist3 = magnitude2D(p1 - p3);
		// compute circumradius (proxy of the max distance Fermat-Torricelli point in the triangle)
		// Simple meaning = if you pick any random point in the domain, if it is less than a circumradius away in aerial view
		// from any vertex of the TIN, then it is very probably inside the TIN, or at least very close from it.
		// From my research it is the quickest way to make this test while keeping (very) good precision, without relying on
		// background grid / octree / subdivision.
		if (isTriangleValid(dist1, dist2, dist3)) {
			double circum = calculateCircumradius(dist1, dist2, dist3);
			if (circum > largest_found && circum<dist1 && circum<dist2 && circum<dist3) largest_found = circum;
		}
	}
	largest_edge_length_ = largest_found;
}

inline int Surface::get_nb_trgls() {
	return int(trgls_.size());
}

inline int Surface::get_nb_pts() {
	return int(pts_.size());
}

inline Triangle Surface::get_triangle(int i) {
	return trgls_[i];
}

inline Vector3 Surface::nearest(Vector3 p1, Vector3 p2, Vector3 ptest) {
	double dist = squaredmagnitude(p1 - ptest);
	double dist2 = squaredmagnitude(p2 - ptest);
	return (dist<=dist2) ? p1 : p2;
}

inline Vector3 Surface::nearest(Vector3 p1, Vector3 p2, Vector3 p3, Vector3 ptest) {

	double dist = squaredmagnitude(p1 - ptest);
	double dist2 = squaredmagnitude(p2 - ptest);
	double dist3 = squaredmagnitude(p3 - ptest);
	if (dist < dist2) {
		if (dist < dist3) {
			return p1;
		}
		return p3;
	}
	else if (dist2 < dist3) {
		return p2;
	}
	return p3;
}

inline Vector3 Surface::get_boundbox_min() {
	return b_min_;
}

inline Vector3 Surface::get_boundbox_max() {
	return b_max_;
}

inline double Surface::get_largest_edge_length() {
	return largest_edge_length_;
}

//*!
//\brief Computes and returns the center of the triangle.
//*/
//inline Vector3 Triangle::Center() const
//{
//	return (pts[0] + pts[1] + pts[2]) / 3.0f;
//}

class Segment
{
protected:
	Vector3 p1_;
	Vector3 p2_;
public:
	Segment();
	Segment(const Vector3& a, const Vector3& b);
	Vector3 start();
	Vector3 end();
};

/*!
\brief Default Constructor.
*/
inline Segment::Segment()
{
}

/*!
\brief 
*/
inline Segment::Segment(const Vector3& a, const Vector3& b)
{
	p1_ = a;
	p2_ = b;
}

inline Vector3 Segment::start()
{
	return p1_;
}

inline Vector3 Segment::end()
{
	return p2_;
}



class Line
{
protected:
	std::vector<Segment> segs_;
	int nb_unique_nodes_;
	std::vector<Vector3> unique_nodes_;
	//std::vector<std::vector<int>> unique_nodes_idx_;
public:
	Line();
	Line(std::vector<Segment> segs);
	void append(Segment seg);
	int size();
	Segment get_seg(int i);
	std::vector<Vector3> get_unique_nodes();
	int get_nb_unique_nodes();
	int get_nb_segs();
};

inline Line::Line(){}

inline Line::Line(std::vector<Segment> segs) {
	segs_ = segs;

	for (int i = 0; i<int(segs_.size()); i++) {

		//int idx = std::find(unique_nodes_.begin(), unique_nodes_.end(), segs_[i].start()) - unique_nodes_.end();

		if (std::find(unique_nodes_.begin(), unique_nodes_.end(), segs_[i].start()) == unique_nodes_.end()) { // node is unique
			unique_nodes_.push_back(segs_[i].start());
		}
		if (std::find(unique_nodes_.begin(), unique_nodes_.end(), segs_[i].end()) == unique_nodes_.end()) { // node is unique
			unique_nodes_.push_back(segs_[i].end());
		}
	}
	nb_unique_nodes_ = int(unique_nodes_.size());
}

inline void Line::append(Segment seg) {
	segs_.push_back(seg);
}

inline int Line::size() {
	return int(segs_.size());
}

inline Segment Line::get_seg(int i) {
	return segs_[i];
}

inline std::vector<Vector3> Line::get_unique_nodes() {
	return unique_nodes_;
}

inline int Line::get_nb_unique_nodes() {
	return nb_unique_nodes_;
}

inline int Line::get_nb_segs() {
	return int(segs_.size());
}

// AABB 3D
class Box
{
protected:
	Vector3 basis;
	Vector3 end;
	Vector3 u, v, w;
	int nu, nv, nw;
	std::vector<double> mat;
	std::vector<double> mat_t;

public:
	inline Box() { }
	explicit Box(const Vector3& basis, const Vector3& u, const Vector3& v, const Vector3& w, const int& nu, const int& nv, const int& nw);


	bool contains(const Vector3&) const;
	Vector3 Center() const;
	double Distance(const Vector3& p) const;
	Vector3 Diagonal() const;
	Vector3 RandomInside() const;
	Vector3 Vertex(int) const;
	Vector3 get_basis() const;
	Vector3 get_end() const;
	Vector3& operator[](int i);
	Vector3 operator[](int i) const;
	Vector3 get_u() const;
	Vector3 get_v() const;
	Vector3 get_w() const;
	int get_nu() const;
	int get_nv() const;
	int get_nw() const;
	void xyz2uvw(const Vector3& pt, int &u, int &v, int &w);

	void Poisson(std::vector<Vector3>& samples, double r, int n) const;
	void Poisson_adapted(std::vector<Vector3>& p, double r, int n, std::vector<Sphere>& spheres) const;
};

/*
\brief Constructor
\param A lower left vertex in world coordinates
\param B upper right vertex in world coordinates
*/
inline Box::Box(const Vector3& basis, const Vector3& u, const Vector3& v, const Vector3& w, const int& nu, const int& nv, const int& nw) : basis(basis), u(u),v(v),w(w),nu(nu),nv(nv),nw(nw)
{
	end = basis + u + v + w;
	mat = { u.x,v.x,w.x,u.y,v.y,w.y,u.z,v.z,w.z };
	double a = mat[0], b = mat[1], c = mat[2], d = mat[3], e = mat[4], f = mat[5], g = mat[6], h = mat[7], i = mat[8];
	double det = a * (e*i - f * h) - b * (d*i - f * g) + c * (d*h - e * g);
	mat_t.push_back(e*i - f * h);
	mat_t.push_back(-(b*i - c * h));
	mat_t.push_back(b*f - c * e);
	mat_t.push_back(-(d*i - f * g));
	mat_t.push_back(a*i - c * g);
	mat_t.push_back(-(a*f - c * d));
	mat_t.push_back(d*h - e * g);
	mat_t.push_back(-(a*h - b * g));
	mat_t.push_back(a*e - b * d);
	for (int iter = 0; iter < int(mat_t.size()); iter++) {
		mat_t[iter] /= det;
	}
}

/*
\brief Returns true if p is inside the box, false otherwise.
\param p world point
*/
inline bool Box::contains(const Vector3& p) const
{
	// we don't know for sure if basis x, y and z are smaller than end x, y and z so we operate this simple test which works in either case:
	double dx = (p.x - basis.x)*(p.x - end.x);
	double dy = (p.y - basis.y)*(p.y - end.y);
	double dz = (p.z - basis.z)*(p.z - end.z);
	return ((dx<=0)&&(dy<=0)&&(dz<=0));
}

/*!
\brief Computes and returns the center of the box.
*/
inline Vector3 Box::Center() const
{
	return (basis + end) / 2.0f;
}

/*!
\brief Returns the diagonal/size of the box.
*/
inline Vector3 Box::Diagonal() const
{
	return (end - basis);
}

/*!
DOESNT WORK AS IT IS NOW
\brief Compute the distance between a point and the box.
\param p point
*/
inline double Box::Distance(const Vector3& p) const
{
	double r = 0.0;
	for (int i = 0; i < 3; i++)
	{
		if (p[i] < basis[i])
		{
			double s = p[i] - basis[i];
			r += s * s;
		}
		else if (p[i] > end[i])
		{
			double s = p[i] - end[i];
			r += s * s;
		}
	}
	return r;
}

/*!
\brief Compute a random point inside a box. Note that this
is not a uniform sampling if the box is not a regular box (width = height = length).

In practice, a uniform sampling doesn't give different results and is way slower, so
we avoid this.

\return a random point inside the box.
*/
inline Vector3 Box::RandomInside() const
{
	Vector3 s = end - basis;
	double randw = Random::Uniform(-1.0f * s[0] / 2.0f, s[0] / 2.0f);
	double randh = Random::Uniform(-1.0f * s[1] / 2.0f, s[1] / 2.0f);
	double randl = Random::Uniform(-1.0f * s[2] / 2.0f, s[2] / 2.0f);
	return (end + basis) / 2.0f + Vector3(randw, randh, randl);
}

/*
\brief Get one of the vertex of the box.
*/
inline Vector3 Box::Vertex(int i) const
{
	if (i == 0)
		return basis;
	return end;
}

/*
\brief Get bottom left vertex in world coordinates
*/
inline Vector3 Box::get_basis() const
{
	return basis;
}

/*
\brief Get top right vertex in world coordinates
*/
inline Vector3 Box::get_end() const
{
	return end;
}

inline Vector3 Box::get_u() const
{
	return u;
}

inline Vector3 Box::get_v() const
{
	return v;
}

inline Vector3 Box::get_w() const
{
	return w;
}

inline int Box::get_nu() const
{
	return nu;
}

inline int Box::get_nv() const
{
	return nv;
}

inline int Box::get_nw() const
{
	return nw;
}

/*
\brief Access box vertex by reference
*/
inline Vector3& Box::operator[](int i)
{
	if (i == 0)
		return basis;
	return end;
}

/*
\brief Access box vertex by const value
*/
inline Vector3 Box::operator[](int i) const
{
	if (i == 0)
		return basis;
	return end;
}

inline void Box::xyz2uvw(const Vector3& pt, int &up, int &vp, int &wp) {

	// 1: normalize vector

	Vector3 ptn;
	ptn.x = mat_t[0] * (pt.x - basis.x) + mat_t[1] * (pt.y - basis.y) + mat_t[2] * (pt.z - basis.z);
	ptn.y = mat_t[3] * (pt.x - basis.x) + mat_t[4] * (pt.y - basis.y) + mat_t[5] * (pt.z - basis.z);
	ptn.z = mat_t[6] * (pt.x - basis.x) + mat_t[7] * (pt.y - basis.y) + mat_t[8] * (pt.z - basis.z);

	// 2: get u/v/w coordinates in normalized domain

	up = ((int)(ceil(ptn.x*nu))) - 1; // values are ceiled (rounded up)
	vp = ((int)(ceil(ptn.y*nv))) - 1;
	wp = ((int)(ceil(ptn.z*nw))) - 1;
}

/*!
\brief Compute a Poisson sphere distribution inside a box.
This function uses a simple O(n<SUP>3</SUP>) dart throwing algorithm.
\param array of existing samples (possibly empty)
\param r Radius of the sphere.
\param n Number of candidate points.
*/
inline void Box::Poisson(std::vector<Vector3>& p, double r, int n) const
{
	double c = 4.0 * r * r;
	for (int i = 0; i < n; i++)
	{
		Vector3 t = RandomInside();
		bool hit = false;
		for (int j = 0; j < int(p.size()); j++)
		{
			if (squaredmagnitude(t - p.at(j)) < c)
			{
				hit = true;
				break;
			}
		}
		if (! hit)
			p.push_back(t);
	}
}


/*!
\brief Compute a Poisson sphere distribution inside a box.
This function uses a simple O(n<SUP>3</SUP>) dart throwing algorithm.
\param array of existing samples (possibly empty)
\param r Radius of the sphere.
\param n Number of candidate points.
\param spheres list of the sphere in which we don't want to add sampling points
*/
inline void Box::Poisson_adapted(std::vector<Vector3>& p, double r, int n, std::vector<Sphere>& spheres) const
{
	double c = 4.0 * r * r;
	for (int i = 0; i < n; i++)
	{
		Vector3 t = RandomInside();
		bool hit = false;
		for (int j = 0; j < int(spheres.size()); j++) {
			if (spheres[j].Contains(t)) {
				hit = true;
				break;
			}
		}
		for (int j = 0; j < int(p.size()); j++)
		{
			if (squaredmagnitude(t - p.at(j)) < c)
			{
				hit = true;
				break;
			}
		}
		if (! hit)
			p.push_back(t);
	}
}

// AABB 2D
class Box2D
{
protected:
	Vector2 a;
	Vector2 b;

public:
	explicit Box2D();
	explicit Box2D(const Vector2& A, const Vector2& B);
	explicit Box2D(const Vector2& C, double R);
	explicit Box2D(const Box& b);

	bool Contains(const Vector2&) const;
	bool Intersect(const Box2D& box) const;
	double Distance(const Vector2& p) const;
	void Poisson(std::vector<Vector2>& p, double r, int n) const;
	Vector2 RandomInside() const;

	Vector2 Vertex(int i) const;
	Vector2 Center() const;
	Vector2 get_basis() const;
	Vector2 get_end() const;
	//Box ToBox(double zMin, double zMax) const;
	Vector2& operator[](int i);
	Vector2 operator[](int i) const;
};

/*
\brief Default Constructor
*/
inline Box2D::Box2D()
{
	a = Vector2(0);
	b = Vector2(0);
}

/*
\brief Constructor
\param A lower left vertex in world coordinates
\param B upper right vertex in world coordinates
*/
inline Box2D::Box2D(const Vector2& A, const Vector2& B) : a(A), b(B)
{
}

/*
\brief Constructor
\param C box center
\param R radius
*/
inline Box2D::Box2D(const Vector2& C, double R)
{
	Vector2 RR = Vector2(R);
	a = C - RR;
	b = C + RR;
}

/*!
\brief Constructor from a 3D box.
\param box the box
*/
inline Box2D::Box2D(const Box& box)
{
	a = Vector2(box.Vertex(0));
	b = Vector2(box.Vertex(1));
}

/*
\brief Returns true if p is inside the box, false otherwise.
\param p world point
*/
inline bool Box2D::Contains(const Vector2& p) const
{
	return (p > a && p < b);
}

/*!
\brief Check if the 2D box intersects another 2D box.
\param box argument box.
*/
inline bool Box2D::Intersect(const Box2D& box) const
{
	if (((a[0] >= box.b[0]) || (a[1] >= box.b[1]) || (b[0] <= box.a[0]) || (b[1] <= box.a[1])))
		return false;
	else
		return true;
}

/*:
\brief Compute the distance between a point and the box.
\param p point
*/
inline double Box2D::Distance(const Vector2 & p) const
{
	double r = 0.0;
	for (int i = 0; i < 2; i++)
	{
		if (p[i] < a[i])
		{
			double s = p[i] - a[i];
			r += s * s;
		}
		else if (p[i] > b[i])
		{
			double s = p[i] - b[i];
			r += s * s;
		}
	}
	return r;
}

/*
\brief Get one of the vertex of the box.
*/
inline Vector2 Box2D::Vertex(int i) const
{
	if (i == 0)
		return a;
	return b;
}

/*
\brief Compute the box center.
*/
inline Vector2 Box2D::Center() const
{
	return (a + b) / 2;
}

/*
\brief Get bottom left vertex in world coordinates
*/
inline Vector2 Box2D::get_basis() const
{
	return a;
}

/*
\brief Get top right vertex in world coordinates
*/
inline Vector2 Box2D::get_end() const
{
	return b;
}

/*
\brief Transform a Box2 in a Box.
\param yMin altitude of the first vertex for the new Box
\param yMax altitude of the second vertex for the new Box
//*/
//inline Box Box2D::ToBox(double yMin, double yMax) const
//{
//	return Box(a.ToVector3(yMin), b.ToVector3(yMax));
//}

/*
\brief Access box vertex by reference
*/
inline Vector2& Box2D::operator[](int i)
{
	if (i == 0)
		return a;
	return b;
}

/*
\brief Access box vertex by const value
*/
inline Vector2 Box2D::operator[](int i) const
{
	if (i == 0)
		return a;
	return b;
}

/*!
\brief Compute a random point inside a box. Note that this
is not a uniform sampling if the box is not a regular box (width = height = length).

In practice, a uniform sampling doesn't give different results and is way slower, so
we avoid this.

\return a random point inside the box.
*/
inline Vector2 Box2D::RandomInside() const
{
	Vector2 s = b - a;
	double randw = Random::Uniform(-1.0f * s[0] / 2.0f, s[0] / 2.0f);
	double randh = Random::Uniform(-1.0f * s[1] / 2.0f, s[1] / 2.0f);
	return (a + b) / 2.0f + Vector2(randw, randh);
}

/*!
\brief Compute a Poisson sphere distribution inside a box.
This function uses a simple O(n<SUP>3</SUP>) dart throwing algorithm.
\param array of existing samples (possibly empty)
\param r Radius of the sphere.
\param n Number of candidate points.
*/
inline void Box2D::Poisson(std::vector<Vector2>& p, double r, int n) const
{
	double c = 4.0 * r * r;
	for (int i = 0; i < n; i++)
	{
		Vector2 t = RandomInside();
		bool hit = false;
		for (int j = 0; j < int(p.size()); j++)
		{
			if (squaredmagnitude(t - p.at(j)) < c)
			{
				hit = true;
				break;
			}
		}
		if (hit == false)
			p.push_back(t);
	}
}


// Circle2. A 2D Circle geometric element.
class Circle2
{
protected:
	Vector2 center;
	double radius;

public:
	Circle2(const Vector2& c, double r);

	Vector2 RandomOn() const;
	Vector2 Center() const;
	double Radius() const;
	bool Contains(const Vector2& p) const;
};

/*!
\brief Constructor
\param c center
\param r radius
*/
inline Circle2::Circle2(const Vector2& c, double r)
{
	center = c;
	radius = r;
}

/*!
\brief Compute a random point on the circle, uniformly.
*/
inline Vector2 Circle2::RandomOn() const
{
	double u = Random::Uniform(-radius, radius);
	double v = Random::Uniform(-radius, radius);
	double s = u * u + v * v;

	double rx = (u * u - v * v) / s;
	double ry = 2.0f * u * v / s;
	return center + Vector2(rx, ry) * radius;
}

/*!
\brief Returns the circle center.
*/
inline Vector2 Circle2::Center() const
{
	return center;
}

/*!
\brief Returns the circle radius.
*/
inline double Circle2::Radius() const
{
	return radius;
}

/*!
\brief Check if a given point lies inside the circle.
*/
inline bool Circle2::Contains(const Vector2& p) const
{
	return (magnitude(p - center) < radius);
}


// Circle. A 3D Circle geometric element.
class Circle
{
protected:
	Vector3 center;
	Vector3 normal;
	double radius;

public:
	Circle(const Vector3& c, const Vector3& n, double r);

	Vector3 Center() const;
	Vector3 Normal() const;
	double Radius() const;
	bool Intersect(const Ray& r, double& t) const;
};

/*
\brief
*/
inline Circle::Circle(const Vector3& c, const Vector3& n, double r)
{
	center = c;
	normal = n;
	radius = r;
}

/*
\brief
*/
inline Vector3 Circle::Center() const
{
	return center;
}

/*
\brief
*/
inline Vector3 Circle::Normal() const
{
	return normal;
}

/*
\brief
*/
inline double Circle::Radius() const
{
	return radius;
}

/*!
\brief Compute the ray-disc intersection.
\param ray The ray.
\param t Intersection depth.
*/
inline bool Circle::Intersect(const Ray& ray, double& t) const
{
	double e = Dot(normal, ray.d);
	if (fabs(e) < 1e-6f)
		return false;
	t = Dot(center - ray.o, normal) / e;
	if (t < 0.0f)
		return false;
	Vector3 p = ray(t) - center;
	if (Dot(p, p) > radius * radius)
		return false;
	return true;
}


/*!
\brief Compute a random point inside a sphere, uniformly.
*/
inline Vector3 Sphere::RandomInside() const
{
	return center + Vector3(Random::Uniform(-radius, radius), Random::Uniform(-radius, radius), Random::Uniform(-radius, radius));
}

/*!
\brief
*/
inline Vector3 Sphere::RandomSurface() const
{
	return Vector3(0);
}

/*!
\brief Returns the sphere center.
*/
inline Vector3 Sphere::Center() const
{
	return center;
}

/*!
\brief Returns the sphere radius.
*/
inline double Sphere::Radius() const
{
	return radius;
}

/*!
\brief Compute a Poisson sphere distribution inside a box.
This function uses a simple O(n<SUP>3</SUP>) dart throwing algorithm.
\param array of existing samples (possibly empty)
\param r Radius of the sphere.
\param n Number of candidate points.
*/
inline void Sphere::Poisson(std::vector<Vector3>& p, double r, int n) const
{
	double c = 4.0 * r * r;
	for (int i = 0; i < n; i++)
	{
		Vector3 t = RandomInside();
		bool hit = false;
		for (int j = 0; j < int(p.size()); j++)
		{
			if (squaredmagnitude(t - p.at(j)) < c)
			{
				hit = true;
				break;
			}
		}
		if (hit == false)
			p.push_back(t);
	}
}

/*!
\brief Compute a Poisson sphere distribution inside a box.
This function uses a simple O(n<SUP>3</SUP>) dart throwing algorithm.
\param array of existing samples (possibly empty)
\param r Radius of the sphere.
\param n Number of candidate points.
*/
inline void Sphere::Poisson_adapted(std::vector<Vector3>& p, double r, int n, std::vector<Sphere>& spheres) const
{
	double c = 4.0 * r * r;
	for (int i = 0; i < n; i++)
	{
		Vector3 t = RandomInside();
		bool hit = false;
		for (int j = 0; j < int(spheres.size()); j++) {
			if (spheres[j].Contains(t)) {
				hit = true;
				break;
			}
		}
		for (int j = 0; j < int(p.size()); j++)
		{
			if (squaredmagnitude(t - p.at(j)) < c)
			{
				hit = true;
				break;
			}
		}
		if (hit == false)
			p.push_back(t);
	}
}


// ScalarField2D. Represents a 2D field (nx * ny) of scalar values bounded in world space. Can represent a heightfield.
class ScalarField2D
{
protected:
	Box2D box;
	int nx, ny;
	std::vector<double> values;

public:
	/*
	\brief Default Constructor
	*/
	inline ScalarField2D() : nx(0), ny(0)
	{
		// Empty
	}

	/*
	\brief Constructor
	\param nx size in x axis
	\param ny size in z axis
	\param bbox bounding box of the domain in world coordinates
	*/
	inline ScalarField2D(int nx, int ny, const Box2D& bbox) : box(bbox), nx(nx), ny(ny)
	{
		values.resize(size_t(nx * ny));
	}

	/*
	\brief Constructor
	\param nx size in x axis
	\param ny size in y axis
	\param bbox bounding box of the domain
	\param value default value of the field
	*/
	inline ScalarField2D(int nx, int ny, const Box2D& bbox, double value) : box(bbox), nx(nx), ny(ny)
	{
		values.resize(nx * ny);
		Fill(value);
	}

	/*!
	\brief Constructor
	\param nx size in x axis
	\param ny size in y axis
	\param bbox bounding box of the domain
	\param value all values for the field
	*/
	inline ScalarField2D(int nx, int ny, const Box2D& bbox, const std::vector<double>& vals) : ScalarField2D(nx, ny, bbox)
	{
		for (int i = 0; i < int(vals.size()); i++)
			values[i] = vals[i];
	}

	/*
	\brief copy constructor
	\param field Scalarfield2D to copy
	*/
	inline ScalarField2D(const ScalarField2D& field) : ScalarField2D(field.nx, field.ny, field.box)
	{
		for (unsigned int i = 0; i < int(values.size()); i++)
			values[i] = field.values[i];
	}

	/*
	\brief Destructor
	*/
	inline ~ScalarField2D()
	{
	}

	/*
	\brief Compute the gradient for the vertex (i, j)
	*/
	inline Vector2 Gradient(int i, int j) const
	{
		Vector2 ret;
		double cellSizeX = (box.Vertex(1).x - box.Vertex(0).x) / (nx - 1);
		double cellSizeY = (box.Vertex(1).y - box.Vertex(0).y) / (ny - 1);

		// X Gradient
		if (i == 0)
			ret.x = (Get(i + 1, j) - Get(i, j)) / cellSizeX;
		else if (i == ny - 1)
			ret.x = (Get(i, j) - Get(i - 1, j)) / cellSizeX;
		else
			ret.x = (Get(i + 1, j) - Get(i - 1, j)) / (2.0f * cellSizeX);

		// Y Gradient
		if (j == 0)
			ret.y = (Get(i, j + 1) - Get(i, j)) / cellSizeY;
		else if (j == nx - 1)
			ret.y = (Get(i, j) - Get(i, j - 1)) / cellSizeY;
		else
			ret.y = (Get(i, j + 1) - Get(i, j - 1)) / (2.0f * cellSizeY);

		return ret;
	}

	/*
	\brief Normalize this field
	*/
	inline void NormalizeField()
	{
		double min = Min();
		double max = Max();
		for (int i = 0; i < ny * nx; i++)
			values[i] = (values[i] - min) / (max - min);
	}

	/*
	\brief Return the normalized version of this field
	*/
	inline ScalarField2D Normalized() const
	{
		ScalarField2D ret(*this);
		double min = Min();
		double max = Max();
		for (int i = 0; i < ny * nx; i++)
			ret.values[i] = (ret.values[i] - min) / (max - min);
		return ret;
	}

	/*!
	\brief Computes and returns the square root of the ScalarField.
	*/
	inline ScalarField2D Sqrt() const
	{
		ScalarField2D ret(*this);
		for (int i = 0; i < int(values.size()); i++)
			ret.values[i] = sqrt(ret.values[i]);
		return ret;
	}

	/*
	\brief Compute a vertex world position including his height.
	*/
	inline Vector3 Vertex(int i, int j) const
	{
		double x = box.Vertex(0).x + i * (box.Vertex(1).x - box.Vertex(0).x) / (nx - 1);
		double y = Get(i, j);
		double z = box.Vertex(0).y + j * (box.Vertex(1).y - box.Vertex(0).y) / (ny - 1);
		return Vector3(z, y, x);
	}

	/*
	\brief Compute a vertex world position including his height.
	*/
	inline Vector3 Vertex(const Vector2i& v) const
	{
		double x = box.Vertex(0).x + v.x * (box.Vertex(1).x - box.Vertex(0).x) / (nx - 1);
		double y = Get(v.x, v.y);
		double z = box.Vertex(0).y + v.y * (box.Vertex(1).y - box.Vertex(0).y) / (ny - 1);
		return Vector3(z, y, x);
	}

	/*
	\brief Get Vertex world position by performing bilinear interpolation.
	\param v world position in 2D
	*/
	inline Vector3 Vertex(const Vector2& v) const
	{
		return Vector3(v.x, GetValueBilinear(v), v.y);
	}

	/*!
	\brief Check if a point lies inside the bounding box of the field.
	*/
	inline bool Inside(const Vector2& p) const
	{
		Vector2 q = p - box.Vertex(0);
		Vector2 d = box.Vertex(1) - box.Vertex(0);

		double u = q[0] / d[0];
		double v = q[1] / d[1];

		int j = int(u * (nx - 1));
		int i = int(v * (ny - 1));

		return Inside(i, j);
	}

	/*!
	\brief Check if a point lies inside the bounding box of the field.
	*/
	inline bool Inside(int i, int j) const
	{
		if (i < 0 || i >= nx || j < 0 || j >= ny)
			return false;
		return true;
	}

	/*!
	\brief Check if a point lies inside the bounding box of the field.
	*/
	inline bool Inside(const Vector2i& v) const
	{
		if (v.x < 0 || v.x >= nx || v.y < 0 || v.y >= ny)
			return false;
		return true;
	}

	/*!
	\brief Utility.
	*/
	inline Vector2i ToIndex2D(const Vector2& p) const
	{
		Vector2 q = p - box.Vertex(0);
		Vector2 d = box.Vertex(1) - box.Vertex(0);

		double u = q[0] / d[0];
		double v = q[1] / d[1];

		int j = int(u * (nx - 1));
		int i = int(v * (ny - 1));

		return Vector2i(i, j);
	}

	/*!
	\brief Utility.
	*/
	inline void ToIndex2D(int index, int& i, int& j) const
	{
		i = index / nx;
		j = index % nx;
	}

	/*!
	\brief Utility.
	*/
	inline Vector2i ToIndex2D(int index) const
	{
		return Vector2i(index / nx, index % nx);
	}

	/*!
	\brief Utility.
	*/
	inline int ToIndex1D(const Vector2i& v) const
	{
		return v.x * nx + v.y;
	}

	/*!
	\brief Utility.
	*/
	inline int ToIndex1D(int i, int j) const
	{
		return i * nx + j;
	}

	/*!
	\brief Returns the value of the field at a given coordinate.
	*/
	inline double Get(int row, int column) const
	{
		int index = ToIndex1D(row, column);
		return values[index];
	}

	/*!
	\brief Returns the value of the field at a given coordinate.
	*/
	inline double Get(int index) const
	{
		return values[index];
	}

	/*!
	\brief Returns the value of the field at a given coordinate.
	*/
	inline double Get(const Vector2i& v) const
	{
		int index = ToIndex1D(v);
		return values[index];
	}

	/*!
	\brief Compute the bilinear interpolation at a given world point.
	\param p world point.
	*/
	inline double GetValueBilinear(const Vector2& p) const
	{
		Vector2 q = p - box.Vertex(0);
		Vector2 d = box.Vertex(1) - box.Vertex(0);

		double texelX = 1.0f / double(nx - 1);
		double texelY = 1.0f / double(ny - 1);

		double u = q[0] / d[0];
		double v = q[1] / d[1];

		int i = int(v * (ny - 1));
		int j = int(u * (nx - 1));

		if (!Inside(i, j) || !Inside(i + 1, j + 1))
			return -1.0f;

		double anchorU = j * texelX;
		double anchorV = i * texelY;

		double localU = (u - anchorU) / texelX;
		double localV = (v - anchorV) / texelY;

		double v1 = Get(i, j);
		double v2 = Get(i + 1, j);
		double v3 = Get(i + 1, j + 1);
		double v4 = Get(i, j + 1);

		return (1 - localU) * (1 - localV) * v1
			+ (1 - localU) * localV * v2
			+ localU * (1 - localV) * v4
			+ localU * localV * v3;
	}

	/*!
	\brief Fill all the field with a given value.
	*/
	inline void Fill(double v)
	{
		std::fill(values.begin(), values.end(), v);
	}

	/*!
	\brief Set a given value at a given coordinate.
	*/
	inline void Set(int row, int column, double v)
	{
		values[ToIndex1D(row, column)] = v;
	}

	/*!
	\brief Set a given value at a given coordinate.
	*/
	inline void Set(const Vector2i& coord, double v)
	{
		values[ToIndex1D(coord)] = v;
	}

	/*!
	\brief Set a given value at a given coordinate.
	*/
	inline void Set(int index, double v)
	{
		values[index] = v;
	}

	/*!
	\brief Compute the maximum of the field.
	*/
	inline double Max() const
	{
		if (int(values.size()) == 0)
			return 0.0f;
		double max = values[0];
		for (int i = 1; i < int(values.size()); i++)
		{
			if (values[i] > max)
				max = values[i];
		}
		return max;
	}

	/*!
	\brief Compute the minimum of the field.
	*/
	inline double Min() const
	{
		if (int(values.size()) == 0)
			return 0.0f;
		double min = values[0];
		for (int i = 1; i < int(values.size()); i++)
		{
			if (values[i] < min)
				min = values[i];
		}
		return min;
	}

	/*!
	\brief Compute the average value of the scalarfield.
	*/
	inline double Average() const
	{
		double sum = 0.0f;
		for (int i = 0; i < int(values.size()); i++)
			sum += values[i];
		return sum / int(values.size());
	}

	/*!
	\brief Returns the size of the array.
	*/
	inline Vector2i Size() const
	{
		return Vector2i(nx, ny);
	}

	/*!
	\brief Returns the size of x-axis of the array.
	*/
	inline int SizeX() const
	{
		return nx;
	}

	/*!
	\brief Returns the size of y-axis of the array.
	*/
	inline int SizeY() const
	{
		return ny;
	}

	/*!
	\brief Returns the bottom left corner of the bounding box.
	*/
	inline Vector2 get_basis() const
	{
		return box.Vertex(0);
	}

	/*!
	\brief Returns the top right corner of the bounding box.
	*/
	inline Vector2 get_end() const
	{
		return box.Vertex(1);
	}

	/*!
	\brief Returns the bounding box of the field.
	*/
	inline Box2D GetBox() const
	{
		return box;
	}

	/*!
	\brief Compute the memory used by the field.
	*/
	inline int Memory() const
	{
		return sizeof(ScalarField2D) + sizeof(double) * int(int(values.size()));
	}
};

}
