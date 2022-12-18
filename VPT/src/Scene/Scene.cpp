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
	AddMaterial({ 0.803922f, 0.803922f, 0.803922f }, 0.04f, 1, 0.f, false, false, 1.5f);
	//materials.back().textureIndex = 1;
	AddMaterial({ 0.156863f, 0.803922f, 0.172549f }, 0.04f, 1, 0.f, false, false, 1.5f);
	AddMaterial({ 0.803922f, 0.152941f, 0.152941f }, 0.04f, 1, 0.f, false, false, 1.5f);
	AddMaterial({ 1.f, 1.f, 1.f }, 0.2f, 500, 0.f, false, false, 1.5f);
	materials.back().textureIndex = 1;
	AddMaterial({ 1.f, 0.8f, 0.8f }, 0.2f, 500, 0.f, false, true, 1.33f);
	materials.back().textureIndex = 2;
	AddMaterial({ 1.f, 1.f, 1.f }, 0.2f, 500, 0.f, false, false, 1.46f);
	AddMaterial({ 0.803922f, 0.803922f, 0.803922f }, 0.04f, 1, 0.f, true, false, 1.5f);
	addSphere(glm::vec3(0.f, -20.5f, 0.f), 20.f, 0);
	addSphere(glm::vec3(1.f, 0.0f, 0.f), 0.5f, 2);
	//addSphere(glm::vec3(-1.f, 0.0f, 0.f), 0.5f, 1);
	//
	//glm::mat4 transform = glm::rotate(glm::mat4(1.f), glm::pi<float>() / 3, { 0.f, 1.f, 0.f });
	glm::mat4 transform = glm::translate(glm::mat4(1.f), glm::vec3(-1.f, 0.0f, 0.f));
	transform = glm::scale(glm::mat4(1.f), glm::vec3(0.5f));
	addModel("assets/cube.obj", transform, 3);
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
bool Scene::Intersect(Ray& ray, Intersection& isect) const
{
	bool hit = false;
	for (uint32_t i = 0; i < shapes.size(); i++)
	{
		switch (shapes[i].type)
		{
		case ShapeType::Triangle:
			if (shapes[i].triangle.Intersect(ray, isect))
			{
				hit = true;
				isect.objIdx = i;
			}
			break;
		case ShapeType::Sphere:
			if (shapes[i].sphere.Intersect(ray, isect))
			{
				hit = true;
				isect.objIdx = i;
			}
			break;
		case ShapeType::Plane:
			if (shapes[i].plane.Intersect(ray, isect))
			{
				hit = true;
				isect.objIdx = i;
			}
			break;
		}
	}
	return hit;
}
SurfaceInteraction Scene::getSurfaceProperties(const Ray& ray, const Intersection& isect) const
{
	switch (shapes[isect.objIdx].type)
	{
	case ShapeType::Triangle:
		return shapes[isect.objIdx].triangle.getSurfaceProperties(ray, isect);
		break;
	case ShapeType::Sphere:
		return shapes[isect.objIdx].sphere.getSurfaceProperties(ray, isect);
		break;
	case ShapeType::Plane:
		return shapes[isect.objIdx].plane.getSurfaceProperties(ray, isect);
		break;
	}
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
void Scene::addModel(const std::string& filepath, Transform transform, uint32_t material)
{
	models.emplace_back(new Model(filepath, transform, material));
	std::vector<Triangle> tris = models.back()->GetTriangles();
	for (auto tri : tris)
	{
		shapes.push_back(Shape(ShapeType::Triangle, tri));
	}

}
void Scene::addSphere(glm::vec3 pos, float r, uint32_t material)
{
	shapes.push_back(Shape(ShapeType::Sphere, Sphere(pos, r, material)));
}
void Scene::addPlane(glm::vec3 normal, float dist, uint32_t material)
{
	//Shape plane{ ShapeType::Plane };
	//plane.plane = Plane(normal, dist, material);
	shapes.push_back(Shape(ShapeType::Plane, Plane(normal, dist, material)));
}
glm::vec3 Scene::getSkyColor(const Ray& ray) const
{
	//return glm::vec3(0.f);
;	float phi = std::atan2(-ray.direction.z, ray.direction.x)+glm::pi<float>();
	float theta = std::acos(-ray.direction.y);
	float u = phi * glm::one_over_two_pi<float>();
	float v = theta * glm::one_over_pi<float>();
	
	return textures[skyboxIndex]->sampleImageTexture(u, v);
}