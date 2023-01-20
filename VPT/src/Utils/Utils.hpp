#pragma once
#include <glm/glm.hpp>

namespace Utils {
	static uint32_t ConvertToRGBA(const glm::vec4& color)
	{
		uint8_t r = (uint8_t)(color.r * 255.0f);
		uint8_t g = (uint8_t)(color.g * 255.0f);
		uint8_t b = (uint8_t)(color.b * 255.0f);
		uint8_t a = (uint8_t)(color.a * 255.0f);

		uint32_t result = (a << 24) | (b << 16) | (g << 8) | r;
		return result;
	}
}
//uint RandomUInt(uint& seed)
//{
//	seed ^= seed << 13;
//	seed ^= seed >> 17;
//	seed ^= seed << 5;
//	return seed;
//}
//float RandomFloat(uint& seed) { return RandomUInt(seed) * 2.3283064365387e-10f; }
//static unsigned int seed = 0x12345678;
//
//static unsigned int RandomUInt()
//{
//	seed ^= seed << 13;
//	seed ^= seed >> 17;
//	seed ^= seed << 5;
//	return seed;
//}
//static float RandomFloat() { return RandomUInt() * 2.3283064365387e-10f; }