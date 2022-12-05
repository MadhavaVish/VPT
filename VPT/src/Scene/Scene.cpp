#include "Scene.hpp"
#include <iostream>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp>
Scene::Scene()
{
	AddMaterial({ 0.803922f, 0.152941f, 0.152941f }, 0.f, 1, 0.f, false, false, 1.46);
	AddMaterial({ 0.9f, 0.9f, 0.9f }, 0.0f, 500, 0.f, false, false, 1.46);
	//materials.back().radiance = glm::vec3(1.f);
	//AddMaterial({ 0.9f, 0.9f, 0.9f }, false, true, 1.5, 1.f);
	addSphere(glm::vec3(-0.5f, 0.5f, 0.f), 1.f, 0);
	addSphere(glm::vec3(1.f, 0.5f, 0.f), 0.5f, 1);
	//addSphere(glm::vec3(0.f, -20.5f, 0.0f), 20.f, 1);
	////addWalls
	//frontBack.emplace_back(PlaneXY(3.f, glm::vec3(0.f, 0.f, 1.f), 3));
	//frontBack.emplace_back(PlaneXY(-3.f, glm::vec3(0.f, 0.f, -1.f), 3));
	topBottom.emplace_back(PlaneXZ(1.f, glm::vec3(0.f, 1.f, 0.f), 1));
	//topBottom.emplace_back(PlaneXZ(-3.f, glm::vec3(0.f, -1.f, 0.f), 3));
	//leftRight.emplace_back(PlaneYZ(3.f, glm::vec3(1.f, 0.f, 0.f), 1));
	//leftRight.emplace_back(PlaneYZ(-3.f, glm::vec3(-1.f, 0.f, 0.f), 2));
	
	//glm::mat4 transform = glm::translate(glm::mat4(1.f), glm::vec3(-1.f));
	//transform = glm::scale(transform, glm::vec3(0.4f));
	//addModel("C:\\Users\\madha\\OneDrive\\Desktop\\VPT\\cube.obj", Transform(transform), 0);
}

bool Scene::Intersect(Ray& ray, float& tHit, SurfaceInteraction& intersection) const
{
	bool hit = false;
	
	//if (ray.direction.x < 0)
	//{
	//	if (leftRight[0].Intersect(ray, tHit, intersection))
	//		hit = true;
	//}
	//else
	//{
	//	if (leftRight[1].Intersect(ray, tHit, intersection))
	//		hit = true;
	//}
	if (ray.direction.y < 0)
	{
		if (topBottom[0].Intersect(ray, tHit, intersection))
			hit = true;
	}
	//else
	//{
	//	if (topBottom[1].Intersect(ray, tHit, intersection))
	//		hit = true;
	//}
	//if (ray.direction.z < 0)
	//{
	//	if (frontBack[0].Intersect(ray, tHit, intersection))
	//		hit = true;
	//}
	//else
	//{
	//	if (frontBack[1].Intersect(ray, tHit, intersection))
	//		hit = true;
	//}
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

void Scene::AddMaterial(const glm::vec3 albedo, const float& ks, const float& exponent, const float &radiance, const bool &mirror, const bool &glass, const float &ior)
{
	Material mat;
	mat.albedo = albedo;
	mat.ks = ks;
	mat.exponent = exponent;
	mat.radiance = radiance;
	mat.mirror = mirror;
	mat.glass = glass;
	mat.ior = ior;
	materials.push_back(mat);
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
