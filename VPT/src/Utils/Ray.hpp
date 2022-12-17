#pragma once
#include <glm/glm.hpp>
#include "SurfaceInteraction.hpp"

class Ray {

public:
	Ray() = default;
	Ray(const glm::vec3& o, const glm::vec3& dir) :
		origin(o), direction(dir), invDir(1.f/dir) {
		
	};
	glm::vec3 operator()(float t) const
	{
		return origin + direction * t;
	}


public:
	glm::vec3 origin;
	glm::vec3 direction;
	glm::vec3 invDir;

};


