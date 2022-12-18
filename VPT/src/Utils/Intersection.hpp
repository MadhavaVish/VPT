#pragma once
#include <glm/glm.hpp>

struct Intersection
{
	float     t_hit{ std::numeric_limits<float>::infinity() };
	uint32_t  objIdx{ 0 };
	uint32_t  lightIdx{ 0 };
	glm::vec2 barycentric{ 0.f };
};