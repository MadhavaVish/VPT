#pragma once
#include "../Utils/Transform.hpp"
#include "../Utils/Ray.hpp"

class Shape{
public:
    //Shape() : objToWorld(&identity) {}
    //Shape(const glm::mat4* objectToWorld) : objToWorld(objectToWorld) {};
    virtual bool Intersect(const Ray &ray, float* t_hit) const = 0;

public:
    //const glm::mat4* objToWorld;
    const int* pointer;

};