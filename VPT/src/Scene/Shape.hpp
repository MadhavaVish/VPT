#pragma once
#include "../Utils/Transform.hpp"
#include "../Utils/Ray.hpp"

class Shape{
public:

    Shape(const Transform* objectToWorld);
    virtual ~Shape();
    virtual bool Intersect(const Ray &ray, float t_hit) const = 0;

public:
    const Transform* objToWorld;

};