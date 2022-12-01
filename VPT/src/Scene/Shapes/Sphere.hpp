#pragma once
#include "Shape.hpp"

class Sphere : public Shape {
public:
	Sphere(glm::vec3 pos, float radius, const int materialIdx);
	bool Intersect(Ray& ray, float& tHit, SurfaceInteraction &intersection) const override;

	glm::vec3 position;
	float radius, invr;
};