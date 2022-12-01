#include "Plane.hpp"


Plane::Plane(glm::vec3 n, float d, int materialIdx) : Shape(materialIdx), normal(n), distance(d)
{};

bool Plane::Intersect(Ray& ray, float &tHit, SurfaceInteraction &intersection) const
{
	float t = -(glm::dot(ray.origin, this->normal) + this->distance) / (glm::dot(ray.direction, this->normal));
	if (t < tHit && t > 0)
	{
		tHit = t;
		intersection.hit_normal = normal;
		intersection.material = materialIndex;
		return true;
	}
	return false;
};

bool PlaneXY::Intersect(Ray& ray, float& tHit, SurfaceInteraction& intersection) const
{
	float t = -(z + ray.origin.z) / ray.direction.z;
	if (t < tHit && t>0)
	{
		intersection.hit_normal = normal;
		intersection.material = materialIndex;
		tHit = t;
		return true;
	}
	return false;
}

bool PlaneXZ::Intersect(Ray& ray, float& tHit, SurfaceInteraction& intersection) const
{
	float t = -(y + ray.origin.y) / ray.direction.y;
	if (t < tHit && t>0)
	{
		intersection.hit_normal = normal;
		intersection.material = materialIndex;
		tHit = t;
		return true;
	}
	return false;
}

bool PlaneYZ::Intersect(Ray& ray, float& tHit, SurfaceInteraction& intersection) const
{
	float t = -(x + ray.origin.x) / ray.direction.x;
	if (t < tHit && t>0)
	{
		intersection.hit_normal = normal;
		intersection.material = materialIndex;
		tHit = t;
		return true;
	}
	return false;
}