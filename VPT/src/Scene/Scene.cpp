#include "Scene.hpp"
#include <iostream>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp>
Scene::Scene()
{
	AddMaterial({ 1.f, 0.2f, 0.5f }, false, false, 1.3);
	AddMaterial({ 1.f, 0.2f, 0.2f }, false, false, 1.0);
	AddMaterial({ 0.2f, 1.f, 0.2f }, false, false, 1.0);
	AddMaterial({ 1.f, 1.f, 1.f }, false, false, 1.0);
	AddMaterial({ 0.9f, 0.9f, 0.9f }, true, false, 1.3);
	AddMaterial({ 0.955, 0.637, 0.538 }, true, false, 1.5);
	addSphere(glm::vec3(-1.f, 0.f, 0.f), 0.5f, 4);
	//addSphere(glm::vec3(1.f, 0.f, 0.f), 0.5f, 5);
	////addSphere(glm::vec3(0.f, -20.5f, 0.0f), 20.f, 1);
	////addWalls
	//frontBack.emplace_back(PlaneXY(3.f, glm::vec3(0.f, 0.f, 1.f), 3));
	//frontBack.emplace_back(PlaneXY(-3.f, glm::vec3(0.f, 0.f, -1.f), 3));
	//topBottom.emplace_back(PlaneXZ(3.f, glm::vec3(0.f, 1.f, 0.f), 3));
	//topBottom.emplace_back(PlaneXZ(-3.f, glm::vec3(0.f, -1.f, 0.f), 3));
	//leftRight.emplace_back(PlaneYZ(3.f, glm::vec3(1.f, 0.f, 0.f), 1));
	//leftRight.emplace_back(PlaneYZ(-3.f, glm::vec3(-1.f, 0.f, 0.f), 2));
	
	glm::mat4 transform = glm::translate(glm::mat4(1.f), glm::vec3(0.f));
	addModel("C:\\Users\\madha\\OneDrive\\Desktop\\VPT\\triangle.obj", Transform(transform), 0);
}

bool Scene::Intersect(Ray& ray, float& tHit, SurfaceInteraction& intersection) const
{
	bool hit = false;
	
	/*if (ray.direction.x < 0)
	{
		if (leftRight[0].Intersect(ray, tHit, intersection))
			hit = true;
	}
	else
	{
		if (leftRight[1].Intersect(ray, tHit, intersection))
			hit = true;
	}
	if (ray.direction.y < 0)
	{
		if (topBottom[0].Intersect(ray, tHit, intersection))
			hit = true;
	}
	else
	{
		if (topBottom[1].Intersect(ray, tHit, intersection))
			hit = true;
	}
	if (ray.direction.z < 0)
	{
		if (frontBack[0].Intersect(ray, tHit, intersection))
			hit = true;
	}
	else
	{
		if (frontBack[1].Intersect(ray, tHit, intersection))
			hit = true;
	}*/
	for (size_t i = 0; i < spheres.size(); i++)
	{
		if (spheres[i].Intersect(ray, tHit, intersection))
		{
			hit = true;
		}
	}
	for (size_t i = 0; i < triangles.size(); i++)
	{
		if (triangles[i].Intersect(ray, tHit, intersection))
		{
			hit = true;
		}
	}
	return hit;
}

void Scene::AddMaterial(glm::vec3 albedo, bool metallic, bool glass, float ior)
{
	materials.emplace_back(Material(albedo, metallic, glass, ior));
}
void Scene::addSphere(glm::vec3 pos, float r, int material)
{
	spheres.emplace_back(Sphere(pos, r, material));
}
void Scene::addModel(const std::string& filepath, Transform transform, int material)
{
	models.push_back(Model(filepath, transform, material));
	//int numTriangles = models[models.size()].numTriangles();
	//triangles.resize(triangles.size() + numTriangles);
	models.back().GetTriangles(triangles);
	//std::cout << models.back().numTriangles() << std::endl;
}
