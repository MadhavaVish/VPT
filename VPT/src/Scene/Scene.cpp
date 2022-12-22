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
	textures.push_back(new Texture("assets/plant.jpg"));
	textures.push_back(new Texture("assets/rockAlbedo.jpg"));
	AddMaterial({ 0.803922f, 0.803922f, 0.803922f }, 0.04f, 1, 0.f, false, false, 1.5f);
	//materials.back().textureIndex = 1;
	AddMaterial({ 0.156863f, 0.803922f, 0.172549f }, 0.04f, 1, 0.f, false, true, 1.5f);
	AddMaterial({ 0.803922f, 0.152941f, 0.152941f }, 0.04f, 1, 0.f, false, false, 1.5f);
	AddMaterial({ 1.f, 1.f, 1.f }, 0.2f, 500, 0.f, false, false, 1.5f);
	materials.back().textureIndex = 1;
	AddMaterial({ 1.f, 1.f, 1.f }, 0.2f, 500, 0.f, false, false, 1.33f);
	materials.back().textureIndex = 2;
	AddMaterial({ 1.f, 1.f, 1.f }, 0.2f, 500, 0.f, false, false, 1.46f);
	AddMaterial({ 0.803922f, 0.803922f, 0.803922f }, 0.04f, 1, 0.f, true, false, 1.5f);
	addSphere(glm::vec3(0.f, -200.f, 0.f), 200.f, 6);
	//addSphere(glm::vec3(1.f, 0.f , 0.f), 1.f, 4);
	//addSphere(glm::vec3(-1.f, 0.0f, 0.f), 1.f, 3);
	//
	glm::mat4 transform = glm::rotate(glm::mat4(1.f), -glm::pi<float>() / 2, { 0.f, 1.f, 0.f });
	//glm::mat4 transform = glm::translate(glm::mat4(1.f), glm::vec3(0.5f, 0.f, -3.f));
	//transform = glm::rotate(transform, -glm::pi<float>() / 2, { 0.f, 1.f, 0.f });
	//addModel("assets/rock.obj", transform, 4);
	transform = glm::translate(glm::mat4(3.f), glm::vec3(0.f, 0.f, 0.f));
	addModel("assets/bunny2.obj", transform, 4);
	transform = glm::translate(glm::mat4(1.f), glm::vec3(5.0f, 0.f, 2.4f));
	transform = glm::scale(transform, glm::vec3(0.2f));
	//addModel("assets/bunny2.obj", glm::mat4(1.f), 3);
	//addPlane({ 0.f, 0.f, 1.f }, 3.f, 0);
	//addPlane({ 0.f, 0.f, -1.f }, 8.f, 0);
	//addPlane({ 0.f, 1.f, 0.f }, -0.1f, 0);
	//addPlane({ 0.f, -1.f, 0.f }, 3.f, 3);
	//addPlane({ 1.f, 0.f, 0.f }, 3.f, 1);
	//addPlane({ -1.f, 0.f, 0.f }, 3.f, 2);
	//addModel("assets/Cornell/Left.obj", glm::mat4(1.f), 1);
	//addModel("assets/Cornell/Right.obj", glm::mat4(1.f), 2);
	//addModel("assets/Cornell/Top.obj", glm::mat4(1.f), 0);
	transform = glm::translate(glm::mat4(1.f), glm::vec3(0.0f, 2.6f, 0.f));
	//addModel("assets/Cornell/Bottom.obj", transform, 0);
	//addModel("assets/Cornell/Back.obj", glm::mat4(1.f), 0);
	for (unsigned int i = 0; i < shapes.size(); i++)
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
bool Scene::Intersect(Ray& ray, Intersection& isect) const
{
	//bool hit = false;
	//for (uint32_t i = 0; i < shapes.size(); i++)
	//{
	//	switch (shapes[i].type)
	//	{
	//	case ShapeType::Triangle:
	//		if (shapes[i].triangle.Intersect(ray, isect))
	//		{
	//			hit = true;
	//			isect.objIdx = i;
	//		}
	//		break;
	//	case ShapeType::Sphere:
	//		if (shapes[i].sphere.Intersect(ray, isect))
	//		{
	//			hit = true;
	//			isect.objIdx = i;
	//		}
	//		break;
	//	case ShapeType::Plane:
	//		if (shapes[i].plane.Intersect(ray, isect))
	//		{
	//			hit = true;
	//			isect.objIdx = i;
	//		}
	//		break;
	//	}
	//}
	//return hit;
	return IntersectBVH(ray, isect);
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
	default:
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
	for (auto &tri : tris)
	{
		shapes.push_back(Shape(tri));
	}

}
void Scene::addSphere(glm::vec3 pos, float r, uint32_t material)
{
	shapes.push_back(Shape(Sphere(pos, r, material)));
}
void Scene::addPlane(glm::vec3 normal, float dist, uint32_t material)
{
	shapes.push_back(Shape(Plane(normal, dist, material)));
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

void Scene::BuildBVH()
{
	unsigned int N = static_cast<unsigned int>(shapes.size());
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
	Subdivide(rootNodeIdx);
	auto t2 = Clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
	printf("BVH (%i nodes) constructed in %.2fms.\n", nodesUsed, time_span.count() * 1000);
}
void Scene::Flatten(unsigned int QNodeIdx, unsigned int BNodeIdx, bool isRoot)
{
	//BVHNode& node = bvhNode[BNodeIdx];
	if ((BNodeIdx == rootNodeIdx || QNodeIdx == qrootNodeIdx )&& !isRoot)
	{
		return;
	}
	QBVHNode& qNode = qbvhNodes[QNodeIdx];
	for (int i = 0; i < 4; i++)
	{
		qNode.child[i] = 0;
		qNode.count[i] = 0;
	}
	
	BVHNode& node = bvhNode[BNodeIdx];
	std::vector<unsigned int>* childCandidateIndices = new std::vector<unsigned int>;
	unsigned int adoptedBVHChildren[4] = { 0 };
	int numChildren = 0;
	if (bvhNode[node.leftFirst].isLeaf())
	{
		//adopt directly
		glm::vec3 max = bvhNode[node.leftFirst].aabbMax;
		glm::vec3 min = bvhNode[node.leftFirst].aabbMin;
		qNode.bminx4[numChildren] = min.x;
		qNode.bminy4[numChildren] = min.y;
		qNode.bminz4[numChildren] = min.z;
		qNode.bmaxx4[numChildren] = max.x;
		qNode.bmaxy4[numChildren] = max.y;
		qNode.bmaxz4[numChildren] = max.z;
		qNode.child[numChildren] = bvhNode[node.leftFirst].leftFirst;
		qNode.count[numChildren] = bvhNode[node.leftFirst].triCount;
		adoptedBVHChildren[numChildren] = node.leftFirst;
		numChildren++;
	}
	else
	{
		
		unsigned int child1 = bvhNode[node.leftFirst].leftFirst;
		unsigned int child2 = bvhNode[node.leftFirst].leftFirst +1;
		childCandidateIndices->push_back(child1);
		childCandidateIndices->push_back(child2);
	}
	if (bvhNode[node.leftFirst + 1].isLeaf())
	{
		//adopt directly
		BVHNode& adoptee = bvhNode[node.leftFirst + 1];
		glm::vec3 max = adoptee.aabbMax;
		glm::vec3 min = adoptee.aabbMin;
		qNode.bminx4[numChildren] = min.x;
		qNode.bminy4[numChildren] = min.y;
		qNode.bminz4[numChildren] = min.z;
		qNode.bmaxx4[numChildren] = max.x;
		qNode.bmaxy4[numChildren] = max.y;
		qNode.bmaxz4[numChildren] = max.z;
		qNode.child[numChildren] = adoptee.leftFirst;
		qNode.count[numChildren] = adoptee.triCount;
		adoptedBVHChildren[numChildren] = node.leftFirst + 1;
		numChildren++;
	}
	else
	{
		unsigned int child1 = bvhNode[node.leftFirst+1].leftFirst;
		unsigned int child2 = bvhNode[node.leftFirst+1].leftFirst + 1;
		childCandidateIndices->push_back(child1);
		childCandidateIndices->push_back(child2);
	}
	while (1)
	{
		if (childCandidateIndices->size() < 1)
		{
			break;
		}
		unsigned int largest = 0;
		float area = 0;
		for (int i = 0; i < childCandidateIndices->size(); i++)
		{
			float candidateArea = CalculateNodeCost(bvhNode[(*childCandidateIndices)[i]]);
			if (candidateArea > area)
			{
				area = candidateArea;
				largest = i;
			}
		}
		//adopt largest
		// iflargest is leaf else increment qnodesused
		unsigned int adopted =( *childCandidateIndices)[largest];
		BVHNode& adoptee = bvhNode[adopted];
		if (adoptee.isLeaf()) {
			glm::vec3 max = adoptee.aabbMax;
			glm::vec3 min = adoptee.aabbMin;
			qNode.bminx4[numChildren] = min.x;
			qNode.bminy4[numChildren] = min.y;
			qNode.bminz4[numChildren] = min.z;
			qNode.bmaxx4[numChildren] = max.x;
			qNode.bmaxy4[numChildren] = max.y;
			qNode.bmaxz4[numChildren] = max.z;
			qNode.child[numChildren] = adoptee.leftFirst;
			qNode.count[numChildren] = adoptee.triCount;
			adoptedBVHChildren[numChildren] = adopted;
			numChildren++;
		}
		else
		{
			glm::vec3 max = adoptee.aabbMax;
			glm::vec3 min = adoptee.aabbMin;
			qNode.bminx4[numChildren] = min.x;
			qNode.bminy4[numChildren] = min.y;
			qNode.bminz4[numChildren] = min.z;
			qNode.bmaxx4[numChildren] = max.x;
			qNode.bmaxy4[numChildren] = max.y;
			qNode.bmaxz4[numChildren] = max.z;
			qNode.child[numChildren] = qnodesUsed++;
			qNode.count[numChildren] = 0;
			adoptedBVHChildren[numChildren] = adopted;
			numChildren++;
			
		}
		//get children of adopted node
		//adoptedBVHChildren[numChildren - 1] = adopted;
		if (numChildren == 4)
			break;
		//childCandidateIndices.erase(childCandidateIndices.begin() + largest);
		std::swap((*childCandidateIndices)[largest], (*childCandidateIndices)[childCandidateIndices->size() - 1]);
		childCandidateIndices->pop_back();
		if (bvhNode[adopted].isLeaf())
			continue;
		else
		{
			unsigned int child1 = bvhNode[adopted].leftFirst;
			unsigned int child2 = bvhNode[adopted].leftFirst + 1;
			childCandidateIndices->push_back(child1);
			childCandidateIndices->push_back(child2);
		}

	}
	for (int i = 0; i < 4; i++)
	{
		if (qNode.count[i] == 0)
		{
			Flatten(qNode.child[i], adoptedBVHChildren[i], false);
		}
	}
	//Flatten(qNode.child[0], adoptedBVHChildren[0], false);
	//Flatten(qNode.child[1], adoptedBVHChildren[1], false);
	//Flatten(qNode.child[2], adoptedBVHChildren[2], false);
	//Flatten(qNode.child[3], adoptedBVHChildren[3], false);
	delete childCandidateIndices;
	return;
}

void Scene::BuildQBVH()
{
	unsigned int N = static_cast<unsigned int>(shapes.size());
	// create the BVH node pool
	qbvhNodes = (QBVHNode*)_aligned_malloc(sizeof(QBVHNode) * N * 2, 64);
	//QBVHNode& root = qbvhNodes[qrootNodeIdx];
	//BVHNode& b_root = bvhNode[rootNodeIdx];
	Flatten(qrootNodeIdx, rootNodeIdx, true);
	std::cout << qnodesUsed << " number of nodes in qbvh" << std::endl;
}
void Scene::UpdateNodeBounds(unsigned int nodeIdx)
{
	BVHNode& node = bvhNode[nodeIdx];
	AABB nodeBound;
	for (unsigned int first = node.leftFirst, i = 0; i < node.triCount; i++)
	{
		unsigned int leafTriIdx = triIdx[first + i];
		Shape& leafTri = shapes[leafTriIdx];
		switch (leafTri.type)
		{
		case ShapeType::Triangle:
			nodeBound.grow(leafTri.triangle.getBounds());
			break;
		case ShapeType::Sphere:
			nodeBound.grow(leafTri.sphere.getBounds());
			break;
		case ShapeType::Plane:
			nodeBound.grow(leafTri.plane.getBounds());
			break;
		}
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
			Shape& shape = shapes[triIdx[node.leftFirst + i]];
			glm::vec3 cent;
			switch (shape.type)
			{
			case ShapeType::Triangle:
				cent = shape.triangle.centroid;
				break;
			case ShapeType::Sphere:
				cent = shape.sphere.position;
				break;
			case ShapeType::Plane:
				cent = -shape.plane.normal * shapes[triIdx[i]].plane.distance;
				break;
			}
			boundsMin = std::min(boundsMin, cent[a]);
			boundsMax = std::max(boundsMax, cent[a]);
		}
		if (boundsMin == boundsMax) continue;
		// populate the bins
		Bin bin[BINS];
		float scale = BINS / (boundsMax - boundsMin);
		for (unsigned int i = 0; i < node.triCount; i++)
		{
			Shape& shape = shapes[triIdx[node.leftFirst + i]];
			glm::vec3 cent;
			int binIdx;
			switch (shape.type)
			{
			case ShapeType::Triangle:
				cent = shape.triangle.centroid;
				binIdx = std::min(BINS - 1, (int)((cent[a] - boundsMin) * scale));
				bin[binIdx].triCount++;
				bin[binIdx].bounds.grow(shape.triangle.getBounds());
				break;
			case ShapeType::Sphere:
				cent = shape.sphere.position;
				binIdx = std::min(BINS - 1, (int)((cent[a] - boundsMin) * scale));
				bin[binIdx].triCount++;
				bin[binIdx].bounds.grow(shape.sphere.getBounds());
				break;
			case ShapeType::Plane:
				cent = -shape.plane.normal * shapes[triIdx[i]].plane.distance;
				binIdx = std::min(BINS - 1, (int)((cent[a] - boundsMin) * scale));
				bin[binIdx].triCount++;
				bin[binIdx].bounds.grow(shape.plane.getBounds());
				break;
			}
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
void Scene::Subdivide(unsigned int nodeIdx)
{
	// terminate recursion
	BVHNode& node = bvhNode[nodeIdx];
	// determine split axis using SAHw
	int axis;
	float splitPos;
	float splitCost = FindBestSplitPlane(node, axis, splitPos);
	float nosplitCost = CalculateNodeCost(node);
	if (splitCost >= nosplitCost) return;
	// in-place partition
	int i = node.leftFirst;
	int j = i + node.triCount - 1;
	while (i <= j)
	{
		Shape& shape = shapes[triIdx[i]];
		glm::vec3 cent;
		switch (shape.type)
		{
		case ShapeType::Triangle:
			cent = shape.triangle.centroid;
			break;
		case ShapeType::Sphere:
			cent = shape.sphere.position;
			break;
		case ShapeType::Plane:
			cent = -shape.plane.normal * shapes[triIdx[i]].plane.distance;
			break;
		}
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
	Subdivide(leftChildIdx);
	Subdivide(rightChildIdx);
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
inline float IntersectAABB_SSE(const Ray& ray, const __m128& bmin4, const __m128& bmax4, const float& t_hit)
{
	static __m128 mask4 = _mm_cmpeq_ps(_mm_setzero_ps(), _mm_set_ps(1, 0, 0, 0));
	__m128 t1 = _mm_mul_ps(_mm_sub_ps(_mm_and_ps(bmin4, mask4), ray.O4), ray.rD4);
	__m128 t2 = _mm_mul_ps(_mm_sub_ps(_mm_and_ps(bmax4, mask4), ray.O4), ray.rD4);
	__m128 vmax4 = _mm_max_ps(t1, t2), vmin4 = _mm_min_ps(t1, t2);
	float tmax = std::min(vmax4.m128_f32[0], std::min(vmax4.m128_f32[1], vmax4.m128_f32[2]));
	float tmin = std::max(vmin4.m128_f32[0], std::max(vmin4.m128_f32[1], vmin4.m128_f32[2]));
	if (tmax >= tmin && tmin < t_hit && tmax > 0) return tmin; else return std::numeric_limits<float>::infinity();
}
bool Scene::IntersectBVH(Ray& ray, Intersection& isect) const
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
				const Shape& shape = shapes[triIdx[node->leftFirst + i]];
				switch (shape.type)
				{
				case ShapeType::Triangle:
					if (shape.triangle.Intersect(ray, isect))
					{
						hit = true;
						isect.objIdx = triIdx[node->leftFirst + i];
					}
					break;
				case ShapeType::Sphere:
					if (shape.sphere.Intersect(ray, isect.t_hit))
					{
						hit = true;
						isect.objIdx = triIdx[node->leftFirst + i];
					}
					break;
				case ShapeType::Plane:
					if (shape.plane.Intersect(ray, isect.t_hit))
					{
						hit = true;
						isect.objIdx = triIdx[node->leftFirst + i];
					}
					break;
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
		float dist1 = IntersectAABB_SSE(ray, child1->aabbMin4, child1->aabbMax4, isect.t_hit);
		float dist2 = IntersectAABB_SSE(ray, child2->aabbMin4, child2->aabbMax4, isect.t_hit);
		//float dist1 = IntersectAABB(ray, child1->aabbMin, child1->aabbMax, isect.t_hit);
		//float dist2 = IntersectAABB(ray, child2->aabbMin, child2->aabbMax, isect.t_hit);

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

bool Scene::OcclusionBVH(Ray& ray, float distance) const
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
				const Shape& shape = shapes[triIdx[node->leftFirst + i]];
				switch (shape.type)
				{
				case ShapeType::Triangle:
					if (shape.triangle.Intersect(ray, dis))
					{
						if (dis * dis < distance) return true;
					}
					break;
				case ShapeType::Sphere:
					if (shape.sphere.Intersect(ray, dis))
					{
						if (dis * dis < distance) return true;
					}
					break;
				case ShapeType::Plane:
					if (shape.plane.Intersect(ray, dis))
					{
						if (dis * dis < distance) return true;
					}
					break;
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
		float dist1 = IntersectAABB_SSE(ray, child1->aabbMin4, child1->aabbMax4, dis);
		float dist2 = IntersectAABB_SSE(ray, child2->aabbMin4, child2->aabbMax4, dis);
		//float dist1 = IntersectAABB(ray, child1->aabbMin, child1->aabbMax, dis);
		//float dist2 = IntersectAABB(ray, child2->aabbMin, child2->aabbMax, dis);

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