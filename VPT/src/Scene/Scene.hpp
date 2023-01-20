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
public:
	enum class ShapeType
	{
		Triangle,
		Sphere,
		Plane
	};
	struct Shape
	{
		Shape(Triangle tri) : type(ShapeType::Triangle), triangle(tri) {}
		Shape(Sphere sph) : type(ShapeType::Sphere), sphere(sph) {}
		Shape(Plane pla) : type(ShapeType::Plane), plane(pla) {}
		ShapeType type;
		union
		{
			Triangle triangle;
			Sphere sphere;
			Plane plane;
		};
	};
private:
	struct BVHNode
	{
		union { struct { glm::vec3 aabbMin; unsigned int leftFirst; }; __m128 aabbMin4; };
		union { struct { glm::vec3 aabbMax; unsigned int triCount; }; __m128 aabbMax4; };
		bool isLeaf() { return triCount > 0; }
	};
	__declspec(align(128))struct QBVHNode
	{
		QBVHNode() { minx4 = miny4 = minz4 = maxx4 = maxy4 = maxz4 = _mm_set1_ps(0.f); };
		union {
			struct { float bminx4[4]; };
			__m128 minx4;
		};
		union {
			struct { float bminy4[4]; };
			__m128 miny4;
		};
		union {
			struct { float bminz4[4]; };
			__m128 minz4;
		};
		union {
			struct { float bmaxx4[4]; };
			__m128 maxx4;
		};
		union {
			struct { float bmaxy4[4]; };
			__m128 maxy4;
		};
		union {
			struct { float bmaxz4[4]; };
			__m128 maxz4;
		};
		int child[4], count[4];
	};

	struct Bin { AABB bounds; int triCount = 0; };
public:
	Scene();
	~Scene();
	const bool Intersect(Ray& ray, Intersection& isect) const;
	SurfaceInteraction getSurfaceProperties(const Ray& ray, const Intersection& isect) const;
	const glm::vec3 getSkyColor(const Ray& ray) const;

	void AddMaterial(const glm::vec3 albedo, const float& ks, const float& exponent, const float& radiance, const bool& mirror, const bool& glass, const float& ior);
	void addModel(const std::string& filepath, Transform transform, uint32_t material);
	const bool OcclusionBVH(Ray& ray, float distance) const;
private:
	void BuildBVH();
	void BuildQBVH();
	void Flatten(unsigned int QNodeIdx, unsigned int BNodeIdx, bool isRoot, int depth);
	void Subdivide(unsigned int nodeIdx, int depth);
	void UpdateNodeBounds(unsigned int nodeIdx);
	float FindBestSplitPlane(BVHNode& node, int& axis, float& splitPos);
	float CalculateNodeCost(BVHNode& node);
	const bool IntersectBVH(Ray& ray, Intersection& isect) const;
	const void IntersectQBVH(Ray& ray, Intersection& isect, const unsigned int nodeIdx) const;
public:
	std::vector<Model*> models;
	std::vector<Triangle> Triangles;
	std::vector<Material> materials;
	std::vector<Texture*> textures;
private:
	uint32_t skyboxIndex;
	std::vector<unsigned int> triIdx;
	BVHNode* bvhNode = 0;
	QBVHNode* qbvhNodes = 0;
	unsigned int rootNodeIdx = 0, nodesUsed = 2;
	unsigned int qrootNodeIdx = 0, qnodesUsed = 1;
};
