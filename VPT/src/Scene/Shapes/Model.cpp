#include "Model.hpp"
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"
//#define TINYOBJ_LOADER_OPT_IMPLEMENTATION
//#include "tinyobj_loader_opt.h"
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/hash.hpp>
#include <unordered_map>
#include <iostream>
namespace std
{
    template <>
    struct hash<Model::Vertex> {
        size_t operator()(Model::Vertex const& vertex) const {
            size_t seed = 0;
            hashCombine(seed, vertex.position, vertex.normal, vertex.uv);
            return seed;
        }
    };
}
Model::Model(const std::string& filepath, Transform transform, const uint32_t materialIdx) : transformation(transform), material(materialIdx)
{
	tinyobj::attrib_t attrib;
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;
	std::string warn, err;

	if (!tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, filepath.c_str()))
	{
        std::cout << "pain" << std::endl;
		throw std::runtime_error(warn + err);
	}
    std::unordered_map<Vertex, uint32_t> uniqueVertices;
    glm::mat4 inverseTranspose = glm::transpose(transform.worldToObj);
	for (const auto& shape : shapes)
	{
		for (const auto& index : shape.mesh.indices)
		{
            Vertex vertex{};

            if (index.vertex_index >= 0) {
                vertex.position = {
                    attrib.vertices[3 * index.vertex_index + 0],
                    attrib.vertices[3 * index.vertex_index + 1],
                    attrib.vertices[3 * index.vertex_index + 2],
                };
            }

            if (index.normal_index >= 0) {
                vertex.normal = {
                    attrib.normals[3 * index.normal_index + 0],
                    attrib.normals[3 * index.normal_index + 1],
                    attrib.normals[3 * index.normal_index + 2],
                };
            }

            if (index.texcoord_index >= 0) {
                vertex.uv = {
                    attrib.texcoords[2 * index.texcoord_index + 0],
                    attrib.texcoords[2 * index.texcoord_index + 1],
                };
            }
            if (uniqueVertices.count(vertex) == 0) {
                uniqueVertices[vertex] = static_cast<uint32_t>(vertices.size());
                vertices.push_back(glm::vec3(transform.objToWorld * glm::vec4(vertex.position, 1.f)));
                normals.push_back(glm::normalize(glm::vec3(inverseTranspose * glm::vec4(vertex.normal, 0.f))));
                uvs.push_back(vertex.uv);
            }
            indices.push_back(uniqueVertices[vertex]);
		}
	}
    //std::cout << vertices.size() << std::endl;
}

std::vector<Triangle> Model::GetTriangles() const
{
    std::vector<Triangle> tris;
    tris.reserve(numTriangles());
    for (uint32_t i = 0; i < numTriangles(); i++)
    {
        tris.emplace_back(Triangle(this, i));
    }
    return tris;
}

const bool Triangle::Intersect(Ray& ray, Intersection& isect) const
{
    //if (model->vertices.size() == 0) return false;
    glm::vec3 v0 = model->vertices[v_indices[0]], v1= model->vertices[v_indices[1]], v2= model->vertices[v_indices[2]];
    glm::vec3 A = v1 - v0;
    glm::vec3 B = v2 - v0;
    glm::vec3 pvec = cross(ray.direction, B);
    float determinant = glm::dot(A, pvec);

    float invDeterminant = 1 / determinant;

    glm::vec3 tvec = ray.origin - v0;
    float u = glm::dot(tvec, pvec) * invDeterminant;
    if (u < 0 || u  > 1) return false;

    glm::vec3 qvec = cross(tvec, A);
    float v = glm::dot(ray.direction, qvec) * invDeterminant;
    if (v < 0 || u + v > 1) return false;

    float t = glm::dot(B, qvec) * invDeterminant;
    if (t < isect.t_hit && t > 0)
    {
        isect.t_hit = t;
        isect.barycentric = glm::vec2(u, v);
        return true;
    }
    return false;
}
const bool Triangle::Intersect(Ray& ray, float& t_hit) const
{
    if (model->vertices.size() == 0) return false;
    glm::vec3 v0 = model->vertices[v_indices[0]], v1 = model->vertices[v_indices[1]], v2 = model->vertices[v_indices[2]];
    glm::vec3 A = v1 - v0;
    glm::vec3 B = v2 - v0;
    glm::vec3 pvec = cross(ray.direction, B);
    float determinant = glm::dot(A, pvec);

    float invDeterminant = 1 / determinant;

    glm::vec3 tvec = ray.origin - v0;
    float u = glm::dot(tvec, pvec) * invDeterminant;
    if (u < 0 || u  > 1) return false;

    glm::vec3 qvec = cross(tvec, A);
    float v = glm::dot(ray.direction, qvec) * invDeterminant;
    if (v < 0 || u + v > 1) return false;

    float t = glm::dot(B, qvec) * invDeterminant;
    if (t < t_hit && t > 0)
    {
        t_hit = t;
        //isect.barycentric = glm::vec2(u, v);
        return true;

    }
    return false;
}
SurfaceInteraction Triangle::getSurfaceProperties(const Ray& ray, const Intersection& isect) const
{
    SurfaceInteraction interaction;
    interaction.materialIdx = model->material;
    float u = isect.barycentric.x;
    float v = isect.barycentric.y;
    interaction.hit_normal = (1 - u - v) * model->normals[v_indices[0]] + u * model->normals[v_indices[1]] + v * model->normals[v_indices[2]];
    interaction.uv = (1 - u - v) * model->uvs[v_indices[0]] + u * model->uvs[v_indices[1]] + v * model->uvs[v_indices[2]];
    return interaction;
}
AABB Triangle::getBounds()
{
    AABB bound;
    glm::vec3 v0 = model->vertices[v_indices[0]], v1 = model->vertices[v_indices[1]], v2 = model->vertices[v_indices[2]];
    bound.grow(v0);
    bound.grow(v1);
    bound.grow(v2);
    return bound;
}