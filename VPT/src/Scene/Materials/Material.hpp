#pragma once
#include <glm/glm.hpp>

struct Material {
	Material(glm::vec3 albedo = { 1.f, 0.2f, 0.2f }, bool metallic = false, bool glass = false, float ior = 1.3)
		: albedo(albedo), metallic(metallic), glass(glass), ior(ior) {}
	glm::vec3 albedo;
	bool metallic;
	bool glass;
	float ior;
};