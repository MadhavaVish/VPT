#include "Sphere.hpp"

class Sphere : public Shape {
public:
	Sphere(const Transform* t,glm::vec3 p, float r) : 
		Shape(t),
		position(p), radius(r), invr(1 / r) 
	{}
	bool Intersect(const int idx, Ray& ray) const 
	{
		glm::vec3 oc = ray.origin - this->position;
		float b = glm::dot(oc, ray.direction);
		float c = glm::dot(oc, oc) - this->radius * this->radius;
		float t, d = b * b - c;

		if (d <= 0) return false;
		d = sqrt(d), t = -b - d;
		if (t < ray.time && t > 0)
		{
			return true;
		};
		t = d - b;
		if (t < ray.time && t > 0) {
			return true;
		};
		return false;
	}

	glm::vec3 GetNormal(const glm::vec3 I) const
	{
		return(I - this->position) * invr;
	}
	glm::vec3 GetAlbedo(const glm::vec3 I) const
	{
		return glm::vec3(1, 0.2f, 0.2f);
	}

	glm::vec3 position;
	float radius, invr;
	const Transform* objToWorld;

};