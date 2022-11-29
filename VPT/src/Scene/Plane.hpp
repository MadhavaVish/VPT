#pragma once
#include "Shape.hpp"
#include "../Utils/Ray.hpp"


class Plane : public Shape {
public:

	Plane(glm::vec3 n, float d, int materialIdx);
	bool Intersect(const int idx, Ray& ray) const;
	glm::vec3 GetNormal(const glm::vec3 I) const;
	glm::vec3 GetAlbedo(const glm::vec3 I) const;

	glm::vec3 normal;
	float distance;
};