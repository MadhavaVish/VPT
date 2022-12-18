#pragma once
#include <glm/glm.hpp>
#include "../Utils/Ray.hpp"


struct BVHNode
{
	union { struct { glm::vec3 aabbMin; unsigned int leftFirst; }; __m128 aabbMin4; };
	union { struct { glm::vec3 aabbMax; unsigned int triCount; }; __m128 aabbMax4; };
	bool isLeaf() { return triCount > 0; }
};
void BuildBVH()
{
	// create the BVH node pool
	bvhNode = (BVHNode*)_aligned_malloc(sizeof(BVHNode) * N * 2, 64);
	// populate triangle index array
	for (int i = 0; i < N; i++) triIdx[i] = i;
	// calculate triangle centroids for partitioning
	for (int i = 0; i < N; i++)
		tri[i].centroid = (tri[i].vertex0 + tri[i].vertex1 + tri[i].vertex2) * 0.3333f;
	// assign all triangles to root node
	BVHNode& root = bvhNode[rootNodeIdx];
	root.leftFirst = 0, root.triCount = N;
	UpdateNodeBounds(rootNodeIdx);
	// subdivide recursively
	Timer t;
	Subdivide(rootNodeIdx);
	printf("BVH (%i nodes) constructed in %.2fms.\n", nodesUsed, t.elapsed() * 1000);
}