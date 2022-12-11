#include "Scene.hpp"
#include <iostream>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp>
Scene::Scene()
{
	textures.push_back(new Texture("assets/skydome.hdr"));
	skyboxIndex = 0;
	textures.push_back(new Texture("assets/earth.jpg"));
	textures.push_back(new Texture("assets/image.jpg"));
	AddMaterial({ 0.803922f, 0.803922f, 0.803922f }, 0.04f, 1, 0.f, false, false, 1.5);
	//materials.back().textureIndex = 1;
	AddMaterial({ 0.156863f, 0.803922f, 0.172549f }, 0.04f, 1, 0.f, false, false, 1.5);
	AddMaterial({ 0.803922f, 0.152941f, 0.152941f }, 0.04f, 1, 0.f, false, false, 1.5);
	AddMaterial({ 1.f, 1.f, 1.f }, 0.2f, 500, 0.f, false, false, 1.5);
	materials.back().textureIndex = 1;
	AddMaterial({ 1.f, 0.8f, 0.8f }, 0.2f, 500, 0.f, false, true, 1.33);
	materials.back().textureIndex = 2;
	AddMaterial({ 1.f, 1.f, 1.f }, 0.2f, 500, 0.f, false, false, 1.46);
	AddMaterial({ 0.803922f, 0.803922f, 0.803922f }, 0.04f, 1, 0.f, true, false, 1.5);
	//addSphere(glm::vec3(0.f, 0.0f, 0.f), 1.f, 3);
	addSphere(glm::vec3(1.f, 0.0f, 0.f), 1.3f, 4);
	//addSphere(glm::vec3(0.f, 0.0f, -4.f), 1.f, 2);
	////addWalls
	frontBack.emplace_back(PlaneXY(-3.f, glm::vec3(0.f, 0.f, 1.f), 6));
	frontBack.emplace_back(PlaneXY(3.f, glm::vec3(0.f, 0.f, -1.f), 6));
	topBottom.emplace_back(PlaneXZ(-1.f, glm::vec3(0.f, 1.f, 0.f), 0));
	topBottom.emplace_back(PlaneXZ(3.f, glm::vec3(0.f, -1.f, 0.f), 0));
	leftRight.emplace_back(PlaneYZ(-3.f, glm::vec3(1.f, 0.f, 0.f), 1));
	leftRight.emplace_back(PlaneYZ(3.f, glm::vec3(-1.f, 0.f, 0.f), 2));
	
	glm::mat4 transform = glm::rotate(glm::mat4(1.f), glm::pi<float>() / 3, { 0.f, 1.f, 0.f });
	transform = glm::scale(transform, glm::vec3(-1.f));
	transform = glm::translate(transform, glm::vec3(-2.f, 0.0f, 0.f));
	addModel("assets/pyramid.obj", Transform(glm::mat4(1.f)), 3);
	//transform = glm::scale(glm::mat4(1.f), glm::vec3(1.f, 1.26f, 0.1f));
	//transform = glm::rotate(transform, glm::pi<float>() / 2, { 1.f, 0.f, 0.f });
	//transform = glm::translate(glm::mat4(1.f), glm::vec3(2.f, 0.f, 0.f));
	//addModel("assets/cube.obj", Transform(transform), 0);
}
Scene::~Scene()
{
	for (auto model : models)
	{
		delete model;
	}
	for (auto texture : textures)
	{
		delete texture;
	}
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
	//if (ray.direction.y < 0)
	//{
	//	if (topBottom[0].Intersect(ray, tHit, intersection))
	//		hit = true;
	//}
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
	models.emplace_back(new Model(filepath, transform, material));
	models.back()->GetTriangles(triangles);
}

glm::vec3 Scene::getSkyColor(const Ray& ray) const
{
	float phi = std::atan2(-ray.direction.z, ray.direction.x)+glm::pi<float>();
	float theta = std::acos(-ray.direction.y);
	float u = phi * glm::one_over_two_pi<float>();
	float v = theta * glm::one_over_pi<float>();
	
	return textures[skyboxIndex]->sampleImageTexture(u, v);
}
