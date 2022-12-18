#pragma once
#include "../Utils/Ray.hpp"
#include "Shapes/Sphere.hpp"
#include "Shapes/Plane.hpp"
#include "Shapes/Model.hpp"
#include "Materials/Material.hpp"
#include "../Utils/SurfaceInteraction.hpp"
#include "../Utils/Intersection.hpp"
#include "../Utils/Transform.hpp"
#include <vector>
#include "Materials/Texture.hpp"
class Scene {
private:
	enum class ShapeType
	{
		Triangle,
		Sphere,
		Plane
	};
	struct Shape
	{
		Shape(ShapeType t, Triangle tri) : type(t), triangle(tri) {}
		Shape(ShapeType t, Sphere sph) : type(t), sphere(sph) {}
		Shape(ShapeType t, Plane pla) : type(t), plane(pla) {}
		ShapeType type;
		union
		{
			Triangle triangle;
			Sphere sphere;
			Plane plane;
		};
	};

	struct BVHNode
	{
		union { struct { glm::vec3 aabbMin; unsigned int leftFirst; }; __m128 aabbMin4; };
		union { struct { glm::vec3 aabbMax; unsigned int triCount; }; __m128 aabbMax4; };
		bool isLeaf() { return triCount > 0; }
	};
	struct Bin { AABB bounds; int triCount = 0; };
public:
	Scene();
	~Scene();
	bool Intersect(Ray& ray, Intersection& isect) const;
	SurfaceInteraction getSurfaceProperties(const Ray& ray, const Intersection& isect) const;
	glm::vec3 getSkyColor(const Ray& ray) const;

	void AddMaterial(const glm::vec3 albedo, const float& ks, const float& exponent, const float& radiance, const bool& mirror, const bool& glass, const float& ior);
	void addModel(const std::string& filepath, Transform transform, uint32_t material);
	void addSphere(glm::vec3 pos, float r, uint32_t material);
	void addPlane(glm::vec3 normal, float dist, uint32_t material);
private:
	void BuildBVH();
	void Subdivide(unsigned int nodeIdx);
	void UpdateNodeBounds(unsigned int nodeIdx);
	float FindBestSplitPlane(BVHNode& node, int& axis, float& splitPos);
	float CalculateNodeCost(BVHNode& node);
	bool IntersectBVH(Ray& ray, Intersection& isect) const;
public:
	std::vector<Model*> models;
	std::vector<Shape> shapes;
	std::vector<Material> materials;
	std::vector<Texture*> textures;
private:
	uint32_t skyboxIndex;
	std::vector<unsigned int> triIdx;
	BVHNode* bvhNode = 0;
	unsigned int rootNodeIdx = 0, nodesUsed = 2;
};
