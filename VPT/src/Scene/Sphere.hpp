#pragma once
#include "Shape.hpp"
#include "../Utils/Ray.hpp"

class Sphere : public Shape {
public:
	Sphere(glm::vec3 pos, float radius, const int materialIdx);
	Sphere(const glm::mat4 &t, float r, const int materialIdx);
	bool Intersect(const Ray& ray, float* t_hit) const override;
	glm::vec3 GetNormal(const glm::vec3 intersection) const;
	glm::vec3 GetAlbedo(const glm::vec3 intersection) const;

	glm::vec3 position;
	float radius, invr;
};