#pragma once
#include "glm/glm.hpp"
#include "../../Utils/Ray.hpp"
#include "../../Utils/SurfaceInteraction.hpp"
#include "../../Utils/Intersection.hpp"
#include "../../Utils/Transform.hpp"
#include "../AABB.hpp"
class Plane{
public:

	Plane(const glm::vec3 n, const float d, const uint32_t materialIdx);

	bool Intersect(Ray& ray, float& t_hit) const;
	SurfaceInteraction getSurfaceProperties(const Ray& ray, const Intersection& isect) const;
	AABB getBounds();
	glm::vec3 normal;
	float distance;
private:
	uint32_t material;
	Frame local;
};