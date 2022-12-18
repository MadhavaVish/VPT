#include "Sphere.hpp"
#include "glm/gtc/constants.hpp"
Sphere::Sphere(const glm::vec3 pos, const float radius, const uint32_t materialIdx) : position(pos), radius(radius), invr(1 / radius), material(materialIdx)
{
}
glm::vec2 Sphere::getUVCoords(const glm::vec3& point) const
{
	glm::vec3 dir = glm::normalize(point);
	float phi = std::atan2(-dir.z, dir.x) + glm::pi<float>();
	float theta = std::acos(-dir.y);
	float u = phi * glm::one_over_two_pi<float>();
	float v = theta * glm::one_over_pi<float>();

	return glm::vec2(u, v);

}
bool Sphere::Intersect(Ray& ray, Intersection& isect) const
{
	glm::vec3 oc = ray.origin - position;
	float b = glm::dot(oc, ray.direction);
	float c = glm::dot(oc, oc) - radius * radius;
	float t, d = b * b - c;
	
	if (d <= 0) return false;
	d = sqrt(d), t = -b - d;
	if (t < isect.t_hit && t > 0)
	{
		isect.t_hit = t;
		return true;
	};
	t = d - b;
	if (t < isect.t_hit && t > 0) {
		isect.t_hit = t;
		return true;
	};
	return false;
}

SurfaceInteraction Sphere::getSurfaceProperties(const Ray& ray, const Intersection& isect) const
{
	SurfaceInteraction interaction;
	interaction.materialIdx = material;
	interaction.hit_normal = (ray(isect.t_hit) - position) * invr;
	interaction.uv = getUVCoords(interaction.hit_normal);
	return interaction;
}