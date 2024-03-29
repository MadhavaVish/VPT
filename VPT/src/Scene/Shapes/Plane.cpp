#include "Plane.hpp"
#include <glm/gtx/common.hpp>

Plane::Plane(const glm::vec3 n, const float d, const uint32_t materialIdx) : normal(glm::normalize(n)), distance(d), local(Frame(normal)), material(materialIdx)
{};

bool Plane::Intersect(Ray& ray, float& t_hit) const
{
	float t = -(glm::dot(ray.origin, this->normal) + this->distance) / (glm::dot(ray.direction, this->normal));
	if (t < t_hit && t > 0)
	{
		t_hit = t;
		return true;
	}
	return false;
};

SurfaceInteraction Plane::getSurfaceProperties(const Ray& ray, const Intersection& isect) const
{
	SurfaceInteraction interaction;
	interaction.materialIdx = material;
	interaction.hit_normal = normal;
	glm::vec3 rayLocal = local.ToLocal(rayPnt(ray, isect.t_hit));
	interaction.uv.x = glm::fmod(abs(rayLocal.x), 1.f);
	interaction.uv.y = glm::fmod(abs(rayLocal.y), 1.f);
	return interaction;
}
//static AABB minSize;
//minSize.grow(glm::vec3(0.0001));
AABB Plane::getBounds()
{
	AABB bound;
	glm::vec3 corner = local.m_X * 10.f + local.m_Y * 10.f;
	glm::vec3 center = -normal * distance;
	bound.grow(center - corner);
	bound.grow(center + corner);
	//AABB minSize;
	//minSize.grow(center - glm::vec3(- 0.001f));
	//minSize.grow(center + glm::vec3(0.001f));
	//bound.grow(minSize);
	return bound;
}