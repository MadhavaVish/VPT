#include "Sphere.hpp"

Sphere::Sphere(glm::vec3 pos, float r) : position(pos), radius(r), invr(1 / r) 
{
}
Sphere::Sphere(const glm::mat4 &t, float r) : radius(r), invr(1/r)
{
	position = t * glm::vec4(0.f);
}

bool Sphere::Intersect(const Ray& ray, float* t_hit) const
{
	glm::vec3 oc = ray.origin - this->position;
	float b = glm::dot(oc, ray.direction);
	float c = glm::dot(oc, oc) - this->radius * this->radius;
	float t, d = b * b - c;
	
	if (d <= 0) return false;
	d = sqrt(d), t = -b - d;
	if (t < *t_hit && t > 0)
	{
		*t_hit = t;
		return true;
	};
	t = d - b;
	if (t < *t_hit && t > 0) {
		*t_hit = t;
		return true;
	};
	return false;
}
glm::vec3 Sphere::GetNormal(const glm::vec3 intersection) const
{
	return (intersection - position) * invr;
}

glm::vec3 Sphere::GetAlbedo(const glm::vec3 intersection) const
{
	return glm::vec3(1, 0.8f, 0.2f);
}
