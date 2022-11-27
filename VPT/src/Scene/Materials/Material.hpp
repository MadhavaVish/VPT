#pragma once
#include <glm/glm.hpp>

struct Material {
	glm::vec3 albedo;
	bool metallic;
	bool glass;
	float ior;
};