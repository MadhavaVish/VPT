#pragma once
#include "../../Utils/Transform.hpp"
#include "../../Utils/Ray.hpp"
#include "../../Utils/SurfaceInteraction.hpp"
#include "../Materials/Material.hpp"
#include "glm/glm.hpp"
class Shape{
public:
    Shape(const int materialIdx) : materialIndex(materialIdx) {};
    //Shape() : objToWorld(&identity) {}
    //Shape(const glm::mat4* objectToWorld) : objToWorld(objectToWorld) {};
    virtual bool Intersect(Ray& ray, float& tHit, SurfaceInteraction& intersection) const = 0;

public:
    //const glm::mat4* objToWorld;
    int materialIndex;

};