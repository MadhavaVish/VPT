#pragma once
#include "../Utils/Ray.hpp"
#include "Sphere.hpp"

#include <vector>

class Scene {
public:
	Scene();

	int Intersect(Ray& ray, float* closest) const;
	void addSphere(glm::vec3 pos, float r);

public:
	std::vector<Sphere> spheres;


};