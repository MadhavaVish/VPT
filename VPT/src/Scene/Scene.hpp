#pragma once
#include "../Utils/Ray.hpp"
#include "Shapes/Sphere.hpp"
#include "Shapes/Plane.hpp"
#include "Shapes/Model.hpp"
#include "Materials/Material.hpp"
#include "../Utils/SurfaceInteraction.hpp"
#include "../Utils/Transform.hpp"
#include <vector>
#include "Materials/Texture.hpp"
class Scene {
public:
	Scene();
	bool Intersect(Ray& ray, float &tHit, SurfaceInteraction& intersection) const;
	void AddMaterial(const glm::vec3 albedo, const float& ks, const float& exponent, const float& radiance, const bool& mirror, const bool& glass, const float& ior);
	void addSphere(glm::vec3 pos, float r, int material);
	void addModel(const std::string& filepath, Transform transform, int material);
	glm::vec3 getSkyColor(const Ray& ray) const;
public:
	std::vector<Sphere> spheres;
	std::vector<Model> models;
	std::vector<Triangle> triangles;
	std::vector<PlaneXY> frontBack;
	std::vector<PlaneXZ> topBottom;
	std::vector<PlaneYZ> leftRight;
	std::vector<Material> materials;
	std::vector<Texture*> textures;
private:
	uint32_t skyboxIndex;
};