#include "Triangle.hpp"


//WARNING: This does not currently account for whether the points are defined in clockwise or counter clockwise order.
Triangle::Triangle(vec3 p1, vec3 p2, vec3 p3, const int materialIdx) :
	Shape(materialIdx), p(p1), q(p2), r(p3), normal(normalize(cross(p2-p1,p3-p1))) {};

Triangle::Triangle(vec3 p1, vec3 p2, vec3 p3, vec3 n, const int materialIdx) :
	Shape(materialIdx), p(p1), q(p2), r(p3), normal(n) {};

bool Triangle::Intersect(const int idx, Ray& ray) const
{
	//Using Moller Trumbore algorithm
	vec3 A = this->q - this->p;
	vec3 B = this->r - this->p;
	vec3 pvec = cross(ray.direction, B);
	float determinant = dot(A, pvec);
	//This will save multiplications later, thus more performance.
	float invDeterminant = 1 / determinant;

	vec3 tvec = ray.origin - this->p;
	float u = dot(tvec, pvec) * invDeterminant;
	if (u < 0 || u  > 1) return false;

	vec3 qvec = cross(tvec, A);
	float v = dot(ray.direction, qvec) * invDeterminant;
	if (v < 0 || u + v > 1) return false;

	// float t = dot(B, qvec) * invDeterminant; //It's not necessary to calculate t since we are not currently storing the value.

	return true;
};

glm::vec3 Triangle::GetNormal(const vec3 intersection) const
{
	return this->normal;
};

glm::vec3 Triangle::GetAlbedo(const vec3 intersection) const
{
	return vec3(1, 0.8f, 0.2f);
};
