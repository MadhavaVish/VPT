#pragma once
#include "glm/glm.hpp"
#include "../../Utils/Ray.hpp"
#include "../../Utils/SurfaceInteraction.hpp"
#include "../../Utils/Intersection.hpp"

class Sphere{
public:
	Sphere(const glm::vec3 pos,const float radius, const uint32_t materialIdx);

	bool Intersect(Ray& ray, Intersection &isect) const;
	SurfaceInteraction getSurfaceProperties(const Ray& ray, const Intersection& isect) const;
private:
	glm::vec2 getUVCoords(const glm::vec3& point) const;
private:
	glm::vec3 position;
	float radius, invr;
	uint32_t material;
};