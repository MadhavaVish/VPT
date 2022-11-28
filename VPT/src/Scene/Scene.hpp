#pragma once
#include "../Utils/Ray.hpp"
#include "Sphere.hpp"
#include "Materials/Material.hpp"
#include <vector>

class Scene {
public:
	Scene();

	int Intersect(Ray& ray, float &closest) const;
	void AddMaterial(glm::vec3 albedo, bool metallic, bool glass, float ior);
	void addSphere(glm::vec3 pos, float r, int material);

public:
	std::vector<Sphere> spheres;
	std::vector<Material> materials;

};