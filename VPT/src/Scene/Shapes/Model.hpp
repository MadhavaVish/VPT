#pragma once

#include <glm/glm.hpp>
#include <memory>
#include <vector>
#include <string>
#include <functional>
#include "../../Utils/Transform.hpp"
#include "../../Utils/Ray.hpp"
#include "../../Utils/SurfaceInteraction.hpp"
#include "../../Utils/Intersection.hpp"
template <typename T, typename... Rest>
void hashCombine(std::size_t& seed, const T& v, const Rest&... rest) {
	seed ^= std::hash<T>{}(v)+0x9e3779b9 + (seed << 6) + (seed >> 2);
	(hashCombine(seed, rest), ...);
};

class Model
{
public:
	struct Vertex {
		glm::vec3 position{};
		glm::vec3 normal{};
		glm::vec2 uv{};
		bool operator==(const Vertex& other) const {
			return position == other.position && normal == other.normal && uv == other.uv;
		}
	};
	friend class Triangle;

	Model(const std::string& filepath, Transform transform,const uint32_t materialIdx);

	std::vector<Triangle> GetTriangles() const;
	uint32_t numTriangles() const
	{
		return indices.size() / 3;
	}

private:
	uint32_t material;
	bool shadeSmooth = true;
	Transform transformation;
	std::vector<uint32_t>  indices;
	std::vector<glm::vec3> vertices;
	std::vector<glm::vec3> normals;
	std::vector<glm::vec2> uvs;
};

class Triangle
{
public:
	Triangle(const Model* model,const int triangle) : model(model)
	{
		v_indices = &model->indices[3 * triangle];
		glm::vec3 v0 = model->vertices[v_indices[0]], v1 = model->vertices[v_indices[1]], v2 = model->vertices[v_indices[2]];
		centroid = (v0 + v1 + v2)/3.f;
	}

	bool Intersect(Ray& ray, Intersection& isect) const;
	SurfaceInteraction getSurfaceProperties(const Ray& ray, const Intersection& isect) const;
private:
	const Model* model;
	const uint32_t* v_indices;
	glm::vec3 centroid;
};