#pragma once
#include <glm/glm.hpp>
#include "SurfaceInteraction.hpp"

__declspec(align(64))struct Ray
{
	Ray() { O4 = D4 = rD4 = _mm_set1_ps(1); }
	union {
		struct { glm::vec3 origin; float d1; }; __m128 O4;
	};
	union {
		struct { glm::vec3 direction; float d2; }; __m128 D4;
	};
	union {
		struct { glm::vec3 invDir; float d3; }; __m128 rD4;
	};

};
static inline Ray getRay(const glm::vec3& o, const glm::vec3& dir) 
{
	Ray r;
	r.origin = o;
	r.direction = dir;
	r.invDir = 1.f / dir;
	return r;
}
static void setDir(Ray& ray, const glm::vec3& dir)
{
	ray.direction = dir;
	ray.invDir = 1.f / dir;
}

static glm::vec3 rayPnt(const Ray& ray, const float& t)
{
	return ray.origin + ray.direction * t;
}
