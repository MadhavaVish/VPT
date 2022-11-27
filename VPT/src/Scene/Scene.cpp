#include "Scene.hpp"
#include <iostream>
Scene::Scene()
{
	addSphere(glm::vec3(0.f, 0.f, 0.f), 0.5f);
	addSphere(glm::vec3(0.f, -30.5f, 0.0f), 30.f);
}

int Scene::Intersect(Ray& ray, float *closest) const
{
	int idx = -1;
	float t_hit = *closest;
	for (size_t i = 0; i < spheres.size(); i++)
	{
		if (spheres[i].Intersect(ray, &t_hit))
		{
			if (t_hit < *closest)
			{
				idx = i;
				*closest = t_hit;
			}
		}
	}
	return idx;
}

void Scene::addSphere(glm::vec3 pos, float r)
{
	spheres.emplace_back(Sphere(pos, r));
}
