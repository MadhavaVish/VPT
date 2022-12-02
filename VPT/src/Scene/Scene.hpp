#pragma once
#include "../Utils/Ray.hpp"
#include "Shapes/Sphere.hpp"
#include "Shapes/Plane.hpp"
#include "Shapes/Model.hpp"
#include "Materials/Material.hpp"
#include "../Utils/SurfaceInteraction.hpp"
#include "../Utils/Transform.hpp"
#include <vector>

class Scene {
public:
	Scene();

	bool Intersect(Ray& ray, float &tHit, SurfaceInteraction& intersection) const;
	void AddMaterial(glm::vec3 albedo, bool metallic, bool glass, float ior);
	void addSphere(glm::vec3 pos, float r, int material);
	void addModel(const std::string& filepath, Transform transform, int material);

public:
	std::vector<Sphere> spheres;
	std::vector<Model> models;
	std::vector<Triangle> triangles;
	std::vector<PlaneXY> frontBack;
	std::vector<PlaneXZ> topBottom;
	std::vector<PlaneYZ> leftRight;
	std::vector<Material> materials;
};