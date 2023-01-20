#include "Scene.hpp"
#include <iostream>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp>
#include "AABB.hpp"
#include <immintrin.h>
#include <chrono>
typedef std::chrono::high_resolution_clock Clock;
Scene::Scene()
{
	textures.push_back(new Texture("assets/skydome.hdr"));
	skyboxIndex = 0;
	//textures.push_back(new Texture("assets/DeskScene/bread.jpg"));
	//textures.push_back(new Texture("assets/DeskScene/desk.jpg"));
	//textures.push_back(new Texture("assets/DeskScene/lampstand.jpg"));
	//textures.push_back(new Texture("assets/DeskScene/plant.jpg"));
	//textures.push_back(new Texture("assets/DeskScene/pot.jpg"));
	//AddMaterial({ 1.f, 1.f, 1.f }, 0.2f, 500, 0.f, false, false, 1.5f);
	//materials.back().textureIndex = 1;
	//AddMaterial({ 1.f, 1.f, 1.f }, 0.2f, 500, 0.f, false, false, 1.5f);
	//materials.back().textureIndex = 2;
	//AddMaterial({ 1.f, 1.f, 1.f }, 0.2f, 500, 0.f, false, false, 1.5f);
	//materials.back().textureIndex = 3;
	//AddMaterial({ 1.f, 1.f, 1.f }, 0.2f, 500, 0.f, false, false, 1.5f);
	//materials.back().textureIndex = 4;
	//AddMaterial({ 0.4f, 0.4f, 0.4f }, 0.2f, 500, 0.f, false, false, 1.5f);
	//materials.back().textureIndex = 5;
	//AddMaterial({ 1.f, 1.f, 1.f }, 0.2f, 500, 2.f, false, false, 1.5f);
	//materials.back().textureIndex = 3;
	//addModel("assets/DeskScene/bread.obj", glm::mat4(1.f), 0);
	//addModel("assets/DeskScene/table.obj", glm::mat4(1.f), 1);
	//addModel("assets/DeskScene/plant.obj", glm::mat4(1.f), 3);
	//addModel("assets/DeskScene/pot.obj", glm::mat4(1.f), 4);
	//addModel("assets/DeskScene/lamp.obj", glm::mat4(1.f), 5);
	//addModel("assets/DeskScene/lampstand.obj", glm::mat4(1.f), 2);

	AddMaterial({ 0.7, 0.34, 0.21 }, 0.2f, 500, 0.f, false, true, 1.5f);
	addModel("assets/bunny.obj", glm::mat4(1.f), 0);
	AddMaterial({ 0.156863f, 0.803922f, 0.172549f }, 0.2f, 500, 0.f, false, false, 1.5f);
	AddMaterial({ 0.803922f, 0.152941f, 0.152941f }, 0.2f, 500, 0.f, false, false, 1.5f);
	AddMaterial({ 0.803922f, 0.803922f, 0.803922f }, 0.2f, 500, 0.f, false, false, 1.5f);
	addModel("assets/Cornell/Left.obj", glm::mat4(1.f), 1);
	addModel("assets/Cornell/Right.obj", glm::mat4(1.f), 2);
	addModel("assets/Cornell/Top.obj", glm::mat4(1.f), 3);
	addModel("assets/Cornell/Bottom.obj", glm::mat4(1.f), 3);
	addModel("assets/Cornell/Back.obj", glm::mat4(1.f), 3);
	for (unsigned int i = 0; i < Triangles.size(); i++)
	{
		triIdx.push_back(i);
	}
	BuildBVH();
	BuildQBVH();
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
const bool Scene::Intersect(Ray& ray, Intersection& isect) const
{
	float initial = isect.t_hit;
	IntersectQBVH(ray, isect, qrootNodeIdx);
	if (isect.t_hit < initial)
	{
		return true;
	}
	else
	{
		return false;
	}
	//return IntersectBVH(ray, isect);
}
SurfaceInteraction Scene::getSurfaceProperties(const Ray& ray, const Intersection& isect) const
{
	return Triangles[isect.objIdx].getSurfaceProperties(ray, isect);
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
	for (auto &tri : tris)
	{
		Triangles.push_back(tri);
	}

}

const glm::vec3 Scene::getSkyColor(const Ray& ray) const
{
	//return glm::vec3(0.f);
;	float phi = std::atan2(-ray.direction.z, ray.direction.x)+glm::pi<float>();
	float theta = std::acos(-ray.direction.y);
	float u = phi * glm::one_over_two_pi<float>();
	float v = theta * glm::one_over_pi<float>();
	
	return textures[skyboxIndex]->sampleImageTexture(u, v);
}
static int BVHDepth = 0;
void Scene::BuildBVH()
{
	unsigned int N = static_cast<unsigned int>(Triangles.size());
	// create the BVH node pool
	bvhNode = (BVHNode*)_aligned_malloc(sizeof(BVHNode) * N * 2, 64);
	// populate triangle index array
	//for (unsigned int i = 0; i < N; i++) triIdx[i] = i;
	// assign all triangles to root node
	BVHNode& root = bvhNode[rootNodeIdx];
	root.leftFirst = 0, root.triCount = N;
	UpdateNodeBounds(rootNodeIdx);
	// subdivide recursively
	auto t1 = Clock::now();
	Subdivide(rootNodeIdx,0);
	auto t2 = Clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
	printf("BVH (%i nodes) constructed in %.8fms.\n", nodesUsed, time_span.count() * 1000);
	std::cout << "BVH Depth: " << BVHDepth << std::endl;
}

static int QBVHDepth = 0;
void Scene::Flatten(unsigned int QNodeIdx, unsigned int BNodeIdx, bool isRoot, int depth)
{
	//BVHNode& node = bvhNode[BNodeIdx];
	if ((BNodeIdx == rootNodeIdx )&& !isRoot)
	{
		return;
	}
	if (depth > QBVHDepth) QBVHDepth = depth;
	QBVHNode& qNode = qbvhNodes[QNodeIdx];
	for (int i = 0; i < 4; i++)
	{
		qNode.child[i] = 0;
		qNode.count[i] = 0;
	}
	BVHNode& node = bvhNode[BNodeIdx];

	std::vector<unsigned int> childCandidateIndices;
	childCandidateIndices.reserve(4);

	int numChildren = 0;

	if (bvhNode[node.leftFirst].isLeaf())
	{
		childCandidateIndices.emplace_back(node.leftFirst);
	}
	else
	{
		childCandidateIndices.emplace_back(bvhNode[node.leftFirst].leftFirst);
		childCandidateIndices.emplace_back(bvhNode[node.leftFirst].leftFirst + 1);
	}
	if (bvhNode[node.leftFirst + 1].isLeaf())
	{
		childCandidateIndices.emplace_back(node.leftFirst+1);
	}
	else
	{
		childCandidateIndices.emplace_back(bvhNode[node.leftFirst + 1].leftFirst);
		childCandidateIndices.emplace_back(bvhNode[node.leftFirst + 1].leftFirst + 1);
	}
	std::sort(childCandidateIndices.begin(), childCandidateIndices.end(), [&](auto const& e1, auto const& e2) {
		return CalculateNodeCost(bvhNode[e1]) < CalculateNodeCost(bvhNode[e1]);
	});

	while (childCandidateIndices.size() > 0)
	{
		unsigned int adopted = childCandidateIndices.back();
		childCandidateIndices.pop_back();
		BVHNode& adoptee = bvhNode[adopted];

		glm::vec3 min = adoptee.aabbMin;
		qNode.bminx4[numChildren] = min.x;
		qNode.bminy4[numChildren] = min.y;
		qNode.bminz4[numChildren] = min.z;
		glm::vec3 max = adoptee.aabbMax;
		qNode.bmaxx4[numChildren] = max.x;
		qNode.bmaxy4[numChildren] = max.y;
		qNode.bmaxz4[numChildren] = max.z;

		if (adoptee.isLeaf())
		{
			qNode.child[numChildren] = adoptee.leftFirst;
			qNode.count[numChildren] = adoptee.triCount;
		}
		else
		{
			qNode.child[numChildren] = qnodesUsed++;
			qNode.count[numChildren] = 0;
			Flatten(qNode.child[numChildren], adopted, false, depth + 1);
		}
		numChildren++;
	}
	return;
}

void Scene::BuildQBVH()
{
	unsigned int N = static_cast<unsigned int>(Triangles.size());
	// create the BVH node pool
	qbvhNodes = (QBVHNode*)_aligned_malloc(sizeof(QBVHNode) * N * 2, 128);
	//QBVHNode& root = qbvhNodes[qrootNodeIdx];
	//BVHNode& b_root = bvhNode[rootNodeIdx];
	auto t1 = Clock::now();
	Flatten(qrootNodeIdx, rootNodeIdx, true, 0);
	auto t2 = Clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
	printf("QBVH (%i nodes) constructed in %.8fms.\n", qnodesUsed, time_span.count() * 1000);
	std::cout << "QBVH Depth: " << QBVHDepth << std::endl;
}
void Scene::UpdateNodeBounds(unsigned int nodeIdx)
{
	BVHNode& node = bvhNode[nodeIdx];
	AABB nodeBound;
	for (unsigned int first = node.leftFirst, i = 0; i < node.triCount; i++)
	{
		unsigned int leafTriIdx = triIdx[first + i];
		Triangle& leafTri = Triangles[leafTriIdx];
		nodeBound.grow(leafTri.getBounds());
	}
	node.aabbMin = nodeBound.bmin;
	node.aabbMax = nodeBound.bmax;
}

static const int BINS = 8;
float Scene::FindBestSplitPlane(BVHNode& node, int& axis, float& splitPos)
{
	float bestCost = std::numeric_limits<float>::infinity();
	for (int a = 0; a < 3; a++)
	{
		float boundsMin = std::numeric_limits<float>::infinity(), boundsMax = -std::numeric_limits<float>::infinity();
		for (unsigned int i = 0; i < node.triCount; i++)
		{
			Triangle& tri = Triangles[triIdx[node.leftFirst + i]];
			glm::vec3 cent = tri.centroid;
			boundsMin = std::min(boundsMin, cent[a]);
			boundsMax = std::max(boundsMax, cent[a]);
		}
		if (boundsMin == boundsMax) continue;
		// populate the bins
		Bin bin[BINS];
		float scale = BINS / (boundsMax - boundsMin);
		for (unsigned int i = 0; i < node.triCount; i++)
		{
			Triangle& tri = Triangles[triIdx[node.leftFirst + i]];
			glm::vec3 cent = tri.centroid;
			int binIdx = std::min(BINS - 1, (int)((cent[a] - boundsMin) * scale));
			bin[binIdx].triCount++;
			bin[binIdx].bounds.grow(tri.getBounds());
		}
		// gather data for the 7 planes between the 8 bins
		float leftArea[BINS - 1], rightArea[BINS - 1];
		int leftCount[BINS - 1], rightCount[BINS - 1];
		AABB leftBox, rightBox;
		int leftSum = 0, rightSum = 0;
		for (int i = 0; i < BINS - 1; i++)
		{
			leftSum += bin[i].triCount;
			leftCount[i] = leftSum;
			leftBox.grow(bin[i].bounds);
			leftArea[i] = leftBox.area();
			rightSum += bin[BINS - 1 - i].triCount;
			rightCount[BINS - 2 - i] = rightSum;
			rightBox.grow(bin[BINS - 1 - i].bounds);
			rightArea[BINS - 2 - i] = rightBox.area();
		}
		// calculate SAH cost for the 7 planes
		scale = (boundsMax - boundsMin) / BINS;
		for (int i = 0; i < BINS - 1; i++)
		{
			float planeCost = leftCount[i] * leftArea[i] + rightCount[i] * rightArea[i];
			if (planeCost < bestCost)
				axis = a, splitPos = boundsMin + scale * (i + 1), bestCost = planeCost;
		}
	}
	return bestCost;
}

float Scene::CalculateNodeCost(BVHNode& node)
{
	glm::vec3 e = node.aabbMax - node.aabbMin; // extent of the node
	float surfaceArea = e.x * e.y + e.y * e.z + e.z * e.x;
	return node.triCount * surfaceArea;
}

void Scene::Subdivide(unsigned int nodeIdx, int depth)
{
	// terminate recursion
	BVHNode& node = bvhNode[nodeIdx];
	// determine split axis using SAHw
	int axis;
	float splitPos;
	float splitCost = FindBestSplitPlane(node, axis, splitPos);
	float nosplitCost = CalculateNodeCost(node);
	if (splitCost >= nosplitCost) return;
	if (depth > BVHDepth) BVHDepth = depth;
	// in-place partition
	int i = node.leftFirst;
	int j = i + node.triCount - 1;
	while (i <= j)
	{
		Triangle& tri = Triangles[triIdx[i]];
		glm::vec3 cent = tri.centroid;
		
		if (cent[axis] < splitPos)
			i++;
		else
			std::swap(triIdx[i], triIdx[j--]);
	}
	// abort split if one of the sides is empty
	int leftCount = i - node.leftFirst;
	if (leftCount == 0 || leftCount == node.triCount) return;
	// create child nodes
	int leftChildIdx = nodesUsed++;
	int rightChildIdx = nodesUsed++;
	bvhNode[leftChildIdx].leftFirst = node.leftFirst;
	bvhNode[leftChildIdx].triCount = leftCount;
	bvhNode[rightChildIdx].leftFirst = i;
	bvhNode[rightChildIdx].triCount = node.triCount - leftCount;
	node.leftFirst = leftChildIdx;
	node.triCount = 0;
	UpdateNodeBounds(leftChildIdx);
	UpdateNodeBounds(rightChildIdx);
	// recurse
	Subdivide(leftChildIdx, depth+1);
	Subdivide(rightChildIdx, depth + 1);
}
inline float IntersectAABB(const Ray& ray, const glm::vec3& bmin, const glm::vec3& bmax, const float& t_hit)
{
	float tx1 = (bmin.x - ray.origin.x) * ray.invDir.x, tx2 = (bmax.x - ray.origin.x) * ray.invDir.x;
	float tmin = glm::min(tx1, tx2), tmax = glm::max(tx1, tx2);
	float ty1 = (bmin.y - ray.origin.y) * ray.invDir.y, ty2 = (bmax.y - ray.origin.y) * ray.invDir.y;
	tmin = glm::max(tmin, glm::min(ty1, ty2)), tmax = glm::min(tmax, glm::max(ty1, ty2));
	float tz1 = (bmin.z - ray.origin.z) * ray.invDir.z, tz2 = (bmax.z - ray.origin.z) * ray.invDir.z;
	tmin = glm::max(tmin, glm::min(tz1, tz2)), tmax = glm::min(tmax, glm::max(tz1, tz2));
	if (tmax >= tmin && tmin < t_hit && tmax > 0) return tmin; else return std::numeric_limits<float>::infinity();
}

const bool Scene::IntersectBVH(Ray& ray, Intersection& isect) const
{
	BVHNode* node = &bvhNode[rootNodeIdx], * stack[64];
	unsigned int stackPtr = 0;
	bool hit = false; 
	while (1)
	{
		if (node->isLeaf())
		{
			for (unsigned int i = 0; i < node->triCount; i++)
			{
				const Triangle& tri = Triangles[triIdx[node->leftFirst + i]];
				if (tri.Intersect(ray, isect))
				{
					hit = true;
					isect.objIdx = triIdx[node->leftFirst + i];
				}
				continue;
			}
			if (stackPtr == 0)
			{
				break;
			}
			else node = stack[--stackPtr];
			continue;
		}
		BVHNode* child1 = &bvhNode[node->leftFirst];
		BVHNode* child2 = &bvhNode[node->leftFirst + 1];
		float dist1 = IntersectAABB(ray, child1->aabbMin, child1->aabbMax, isect.t_hit);
		float dist2 = IntersectAABB(ray, child2->aabbMin, child2->aabbMax, isect.t_hit);

		if (dist1 > dist2) { std::swap(dist1, dist2); std::swap(child1, child2); }
		if (dist1 == std::numeric_limits<float>::infinity())
		{
			if (stackPtr == 0)
			{
				break; 
			}
			else node = stack[--stackPtr];
		}
		else
		{
			node = child1;
			if (dist2 != std::numeric_limits<float>::infinity()) stack[stackPtr++] = child2;
		}
	}
	return hit;
}

const __m128 IntersectQuadAABB(Ray& ray, __m128 minx4, __m128 miny4, __m128 minz4, __m128 maxx4, __m128 maxy4, __m128 maxz4, const float& t_hit)
{
	__m128 ox = _mm_set_ps1(ray.origin.x), oy = _mm_set_ps1(ray.origin.y), oz = _mm_set_ps1(ray.origin.z);
	__m128 idx = _mm_set_ps1(ray.invDir.x), idy = _mm_set_ps1(ray.invDir.y), idz = _mm_set_ps1(ray.invDir.z);
	__m128 hit = _mm_set_ps1(t_hit);
	__m128 tx1 = _mm_mul_ps(_mm_sub_ps(minx4, ox), idx);
	__m128 tx2 = _mm_mul_ps(_mm_sub_ps(maxx4, ox), idx);
	__m128 tmin = _mm_min_ps(tx1, tx2), tmax = _mm_max_ps(tx1, tx2);
	__m128 ty1 = _mm_mul_ps(_mm_sub_ps(miny4, oy), idy);
	__m128 ty2 = _mm_mul_ps(_mm_sub_ps(maxy4, oy), idy);
	tmin = _mm_max_ps(tmin, _mm_min_ps(ty1, ty2)), tmax = _mm_min_ps(tmax, _mm_max_ps(ty1, ty2));
	__m128 tz1 = _mm_mul_ps(_mm_sub_ps(minz4, oz), idz);
	__m128 tz2 = _mm_mul_ps(_mm_sub_ps(maxz4, oz), idz);
	tmin = _mm_max_ps(tmin, _mm_min_ps(tz1, tz2)), tmax = _mm_min_ps(tmax, _mm_max_ps(tz1, tz2));
	__m128 dists = _mm_set_ps1(std::numeric_limits<float>::infinity());
	__m128 maxge = _mm_cmpgt_ps(tmax, _mm_set_ps1(0.f));
	__m128 lesst = _mm_cmplt_ps(tmin, hit);
	__m128 maxgemin = _mm_cmpge_ps(tmax, tmin);
	__m128 mask = _mm_and_ps(maxge, _mm_and_ps(lesst, maxgemin));
	return _mm_blendv_ps(dists, tmin, mask);
}
const void Scene::IntersectQBVH(Ray& ray, Intersection& isect, const unsigned int nodeIdx) const
{
	QBVHNode* node = &qbvhNodes[nodeIdx];
	alignas(16) float dists[4];
	_mm_store_ps(dists, IntersectQuadAABB(ray, node->minx4, node->miny4, node->minz4, node->maxx4, node->maxy4, node->maxz4, isect.t_hit));
	int distIndices[4]{ 0,1,2,3 };
	if (dists[0] > dists[1]) std::swap(distIndices[0], distIndices[1]);
	if (dists[2] > dists[3]) std::swap(distIndices[2], distIndices[3]);
	if (dists[0] > dists[2]) std::swap(distIndices[0], distIndices[2]);
	if (dists[1] > dists[3]) std::swap(distIndices[1], distIndices[3]);
	if (dists[1] > dists[2]) std::swap(distIndices[1], distIndices[2]);

	for (int j = 0; j < 4; j++)
	{
		if (dists[distIndices[j]] == std::numeric_limits<float>::infinity()) continue;
		if (node->count[distIndices[j]] > 0)
		{
			for (int i = 0; i < node->count[distIndices[j]]; i++)
			{
				const Triangle& tri = Triangles[triIdx[node->child[distIndices[j]] + i]];
				if (tri.Intersect(ray, isect))
				{
					isect.objIdx = triIdx[node->child[distIndices[j]] + i];
				}
				continue;
			}
		}
		else
		{
			IntersectQBVH(ray, isect, node->child[distIndices[j]]);
		}
	}
	return;
}

const bool Scene::OcclusionBVH(Ray& ray, float distance) const
{
	BVHNode* node = &bvhNode[rootNodeIdx], * stack[64];
	unsigned int stackPtr = 0;
	float dis = distance;
	while (1)
	{
		if (node->isLeaf())
		{
			for (unsigned int i = 0; i < node->triCount; i++)
			{
				const Triangle& tri = Triangles[triIdx[node->leftFirst + i]];
				if (tri.Intersect(ray, dis))
				{
					if (dis * dis < distance) return true;
				}
				continue;
			}
			if (stackPtr == 0)
			{
				break;
			}
			else node = stack[--stackPtr];
			continue;
		}
		BVHNode* child1 = &bvhNode[node->leftFirst];
		BVHNode* child2 = &bvhNode[node->leftFirst + 1];
		float dist1 = IntersectAABB(ray, child1->aabbMin, child1->aabbMax, dis);
		float dist2 = IntersectAABB(ray, child2->aabbMin, child2->aabbMax, dis);

		if (dist1 > dist2) { std::swap(dist1, dist2); std::swap(child1, child2); }
		if (dist1 == std::numeric_limits<float>::infinity())
		{
			if (stackPtr == 0)
			{
				break;
			}
			else node = stack[--stackPtr];
		}
		else
		{
			node = child1;
			if (dist2 != std::numeric_limits<float>::infinity()) stack[stackPtr++] = child2;
		}
	}
	return false;
}
