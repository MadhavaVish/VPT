#pragma once
#include "../Scene/Materials/Material.hpp"
#include <glm/glm.hpp>
struct SurfaceInteraction {
	glm::vec3 hit_normal;
	glm::vec2 uv;
	uint32_t  materialIdx = -1;
};

