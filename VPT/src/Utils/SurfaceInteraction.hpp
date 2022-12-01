#pragma once
#include "../Scene/Shape.hpp"
#include "../Scene/Materials/Material.hpp"

struct SurfaceInteraction {
	glm::vec3 hit_normal;
	int material = 0;
	//Light* light = nullptr;
	glm::vec2 uv;
};