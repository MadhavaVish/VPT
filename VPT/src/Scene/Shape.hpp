#pragma once
#include "../Utils/Transform.hpp"
#include "../Utils/Ray.hpp"
#include "../Utils/SurfaceInteraction.hpp"
#include "Materials/Material.hpp"
#include "glm/glm.hpp"
class Shape{
public:
    Shape(const int materialIdx) : materialIndex(materialIdx) {};
    //Shape() : objToWorld(&identity) {}
    //Shape(const glm::mat4* objectToWorld) : objToWorld(objectToWorld) {};
    virtual bool Intersect(const Ray &ray, float* t_hit) const = 0;

public:
    //const glm::mat4* objToWorld;
    const int* pointer;
    int materialIndex;

};