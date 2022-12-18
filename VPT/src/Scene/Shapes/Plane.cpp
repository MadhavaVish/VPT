#include "Plane.hpp"
#include <glm/gtx/common.hpp>

Plane::Plane(const glm::vec3 n, const float d, const uint32_t materialIdx) : normal(n), distance(d), local(Frame(normal)), material(materialIdx)
{};

bool Plane::Intersect(Ray& ray, Intersection& isect) const
{
	float t = -(glm::dot(ray.origin, this->normal) + this->distance) / (glm::dot(ray.direction, this->normal));
	if (t < isect.t_hit && t > 0)
	{
		isect.t_hit = t;
		return true;
	}
	return false;
};

SurfaceInteraction Plane::getSurfaceProperties(const Ray& ray, const Intersection& isect) const
{
	SurfaceInteraction interaction;
	interaction.materialIdx = material;
	interaction.hit_normal = normal;
	glm::vec3 rayLocal = local.ToLocal(ray(isect.t_hit));
	interaction.uv.x = glm::fmod(abs(rayLocal.x), 1.f);
	interaction.uv.y = glm::fmod(abs(rayLocal.y), 1.f);
	return interaction;
}