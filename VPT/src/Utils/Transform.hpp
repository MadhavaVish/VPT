#pragma once
#include <glm/glm.hpp>
//glm::mat4 identity(1.0f);
class Transform {
	Transform() : transform() {};
	Transform(const glm::mat4 &transform) : transform(transform) {};

public:
	const glm::mat4 transform = glm::mat4(1.f);
};