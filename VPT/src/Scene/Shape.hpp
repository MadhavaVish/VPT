#pragma once
#include "../Utils/Transform.hpp"
#include "../Utils/Ray.hpp"
#include "glm/glm.hpp"

class Shape{
public:
    Shape();
    Shape(const Transform* t) : objToWorld(t){};
    virtual ~Shape();
    virtual bool Intersect(const Ray &ray, float t_hit) const = 0;
    virtual glm::vec3 GetNormal(const glm::vec3 I) const = 0;
    virtual glm::vec3 GetAlbedo(const glm::vec3 I) const = 0;

public:
    const Transform* objToWorld;

};