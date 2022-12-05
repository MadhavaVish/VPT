#include "Sphere.hpp"

Sphere::Sphere(glm::vec3 pos, float r, const int materialIdx) : Shape(materialIdx), position(pos), radius(r), invr(1 / r)
{
}
glm::vec2 Sphere::getUVCoords(const glm::vec3& point) const
{
	glm::vec3 dir = glm::normalize(point);
	float phi = std::atan2(-dir.z, dir.x) + glm::pi<float>();
	float theta = std::acos(-dir.y);
	float u = phi * glm::one_over_two_pi<float>();
	float v = 1 - theta * glm::one_over_pi<float>();

	return glm::vec2(u, v);

}
bool Sphere::Intersect(Ray& ray, float& tHit, SurfaceInteraction& intersection) const
{
	glm::vec3 oc = ray.origin - position;
	float b = glm::dot(oc, ray.direction);
	float c = glm::dot(oc, oc) - radius * radius;
	float t, d = b * b - c;
	
	if (d <= 0) return false;
	d = sqrt(d), t = -b - d;
	if (t < tHit && t > 0)
	{
		tHit = t;
		intersection.hit_normal = (ray(t) - position) * invr;
		intersection.material = materialIndex;
		intersection.uv = getUVCoords(intersection.hit_normal);
		return true;
	};
	t = d - b;
	if (t < tHit && t > 0) {
		tHit = t;
		intersection.hit_normal = (ray(t) - position) * invr;
		intersection.material = materialIndex;
		intersection.uv = getUVCoords(intersection.hit_normal);
		return true;
	};
	return false;
}