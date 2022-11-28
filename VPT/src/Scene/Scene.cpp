#include "Scene.hpp"
#include <iostream>
Scene::Scene()
{
	AddMaterial({ 1.f, 0.2f, 0.5f }, false, false, 1.3);
	AddMaterial({ 0.f, 0.4f, 0.8f }, true, false, 1.3);
	AddMaterial({ 0.9999f, 0.9999f,  0.9999f }, false, true, 1.5);
	addSphere(glm::vec3(0.f, 0.f, 0.f), 0.5f, 0);
	addSphere(glm::vec3(0.f, 0.f, 1.f), 0.5f, 2);
	addSphere(glm::vec3(0.f, -20.5f, 0.0f), 20.f, 1);
}

int Scene::Intersect(Ray& ray, float &closest) const
{
	int idx = -1;
	float t_hit = closest;
	for (size_t i = 0; i < spheres.size(); i++)
	{
		if (spheres[i].Intersect(ray, &t_hit))
		{
			if (t_hit < closest)
			{
				idx = i;
				closest = t_hit;
			}
		}
	}
	return idx;
}

void Scene::AddMaterial(glm::vec3 albedo, bool metallic, bool glass, float ior)
{
	materials.emplace_back(Material(albedo, metallic, glass, ior));
}
void Scene::addSphere(glm::vec3 pos, float r, int material)
{
	spheres.emplace_back(Sphere(pos, r, material));
}

