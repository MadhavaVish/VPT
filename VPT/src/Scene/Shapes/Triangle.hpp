#pragma once
#include "Shape.hpp"

using glm::vec3;

class Triangle : public Shape {
public:

	Triangle(vec3 p1, vec3 p2, vec3 p3, const int materialIdx);
	Triangle(vec3 p1, vec3 p2, vec3 p3, vec3 n, const int materialIdx);
	bool Intersect(const int idx, Ray& ray) const;
	vec3 GetNormal(const vec3 I) const;
	vec3 GetAlbedo(const vec3 I) const;

	vec3 p;
	vec3 q; 
	vec3 r;
	vec3 normal;
};