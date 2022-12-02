#pragma once

#include <glm/glm.hpp>
#include <memory>
#include <vector>
#include <string>
#include <functional>
#include "Shape.hpp"
#include "../../Utils/Transform.hpp"
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
	Model(const std::string& filepath, Transform transform,int materialIdx);
	void GetTriangles(std::vector<Triangle> &triangles) const;
	uint32_t numTriangles() const
	{
		return indices.size() / 3;
	}

private:
	int materialIndex;
	bool shadeSmooth = true;
	Transform transformation;
	std::vector<glm::vec3> vertices;
	std::vector<glm::vec3> normals;
	std::vector<glm::vec2> uvs;
	std::vector<uint32_t>  indices;
};

class Triangle
{
public:
	Triangle(const Model* model, int triangle) : model(model)
	{
		v_indices = &model->indices[3 * triangle];
	}

	bool Intersect(Ray& ray, float& tHit, SurfaceInteraction& intersection) const;

private:
	const Model* model;
	const uint32_t* v_indices;
};