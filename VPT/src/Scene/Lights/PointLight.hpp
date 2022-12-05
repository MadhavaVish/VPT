#pragma once
#include <glm/glm.hpp>

struct PointLight
{
	float intensity;
	glm::vec3 color;
	glm::vec3 position;
};