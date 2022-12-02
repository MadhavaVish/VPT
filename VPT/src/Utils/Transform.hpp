#pragma once
#include <glm/glm.hpp>
class Transform {
public:
	Transform() {};
	Transform(glm::mat4 transform) : objToWorld(transform)
	{
		worldToObj = glm::inverse(objToWorld);
	};

public:
	glm::mat4 objToWorld;
	glm::mat4 worldToObj;
};