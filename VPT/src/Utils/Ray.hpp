#pragma once
#include <glm/glm.hpp>
#include "SurfaceInteraction.hpp"

struct Ray {
	Ray() = default;
	Ray(const glm::vec3& o, const glm::vec3& dir) :
		origin(o), direction(dir), invDir(1.f/dir) {};
	Ray(const glm::vec3& dir) :
		origin({1.f}), direction(dir), invDir(1.f / dir) {};
	glm::vec3 operator()(float t) const
	{
		return origin + direction * t;
	}
	glm::vec3 origin;
	glm::vec3 direction;
	glm::vec3 invDir;

};


