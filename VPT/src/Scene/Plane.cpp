#include "Plane.hpp"


Plane::Plane(glm::vec3 n, float d) : normal(n), distance(d)
{};

bool Plane::Intersect(const int idx, Ray& ray) const
{
	float t = -(glm::dot(ray.origin, this->normal) + this->distance) / (glm::dot(ray.direction, this->normal));
	if (t < ray.time && t > 0) return true;
};

glm::vec3 Plane::GetNormal(const glm::vec3 I) const
{
	return this->normal;
};	

glm::vec3 Plane::GetAlbedo(const glm::vec3 I) const
{
	if (this->normal.y == 1)
	{
		// floor albedo: checkerboard
		int ix = (int)(I.x * 2 + 96.01f);
		int iz = (int)(I.z * 2 + 96.01f);
		// add deliberate aliasing to two tile
		if (ix == 98 && iz == 98) ix = (int)(I.x * 32.01f), iz = (int)(I.z * 32.01f);
		if (ix == 94 && iz == 98) ix = (int)(I.x * 64.01f), iz = (int)(I.z * 64.01f);
		return glm::vec3(((ix + iz) & 1) ? 1 : 0.3f);
	}
	else if (this->normal.z == -1)
	{
		//Let's ignore the logo for now.
		//// back wall: logo
		//static Surface logo("assets/logo.png");
		//int ix = (int)((I.x + 4) * (128.0f / 8));
		//int iy = (int)((2 - I.y) * (64.0f / 3));
		//uint p = logo.pixels[(ix & 127) + (iy & 63) * 128];
		//uint3 i3((p >> 16) & 255, (p >> 8) & 255, p & 255);
		//return glm:vec3(i3) * (1.0f / 255.0f);
	}
	return glm::vec3(0.93f);
};