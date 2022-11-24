#pragma once
#include <glm/glm.hpp>

class Ray {

public:
	Ray() : t_max(std::numeric_limits<float>::infinity()), time(0.f) {};
	Ray(const glm::vec3& o, const glm::vec3& dir, float t_max = std::numeric_limits<float>::infinity(), float t = 0.f) :
		origin(o), direction(dir), t_max(t_max), time(t) {

	};
	glm::vec3 operator()(float t) const
	{
		return origin + direction * t;
	}


public:
	glm::vec3 origin;
	glm::vec3 direction;

	float time;
	mutable float t_max;

};