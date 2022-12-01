#include "Sphere.hpp"

Sphere::Sphere(glm::vec3 pos, float r, const int materialIdx) : Shape(materialIdx), position(pos), radius(r), invr(1 / r)
{
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
		return true;
	};
	t = d - b;
	if (t < tHit && t > 0) {
		tHit = t;
		intersection.hit_normal = (ray(t) - position) * invr;
		intersection.material = materialIndex;
		return true;
	};
	return false;
}