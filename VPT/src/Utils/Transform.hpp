#pragma once
#include <glm/glm.hpp>

class Transform {
	Transform() : transform(glm::mat4(1.0f)) {};
	Transform(glm::mat4 transform) : transform(transform) {};

public:
	glm::mat4 transform;
};