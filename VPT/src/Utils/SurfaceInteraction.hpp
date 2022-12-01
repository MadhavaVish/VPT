#pragma once
#include "../Scene/Materials/Material.hpp"

struct SurfaceInteraction {
	glm::vec3 hit_normal;
	int material = 0;
	glm::vec2 uv;
};