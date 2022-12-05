#pragma once

#include <glm/glm.hpp>

struct AreaLight
{
	uint32_t triangleIndex;
	float Area;
	float Intensity;
	float width;
	float height;
	glm::vec3 color;
};
