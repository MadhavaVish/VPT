#pragma once
#include "Shape.hpp"
#include "../Utils/Ray.hpp"

class Sphere : public Shape {
public:
	Sphere(const Transform* t, glm::vec3 p, float r);
	bool Intersect(const int idx, Ray& ray) const;
	glm::vec3 GetNormal(const glm::vec3 I) const;
	glm::vec3 GetAlbedo(const glm::vec3 I) const;

	glm::vec3 position;
	float radius, invr;
	const Transform* objToWorld;
};