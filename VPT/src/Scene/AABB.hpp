#pragma once
#include <glm/glm.hpp>
struct AABB
{
	glm::vec3 bmin{ std::numeric_limits<float>::infinity() }, bmax{ -std::numeric_limits<float>::infinity() };
	void grow(glm::vec3 p) { bmin = glm::min(bmin, p); bmax = glm::max(bmax, p); }
	void grow(AABB& b) { if (b.bmin.x != std::numeric_limits<float>::infinity()) { grow(b.bmin); grow(b.bmax); } }
	float area()
	{
		glm::vec3 e = bmax - bmin; // box extent
		return e.x * e.y + e.y * e.z + e.z * e.x;
	}
};