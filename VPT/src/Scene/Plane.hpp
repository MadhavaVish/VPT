#pragma once
#include "Shape.hpp"
#include "../Utils/Ray.hpp"


class Plane : public Shape
{
public:

	Plane(glm::vec3 n, float d, int materialIdx);
	bool Intersect(Ray& ray, float& tHit, SurfaceInteraction &intersection) const;
private:

	glm::vec3 normal;
	float distance;
};

class PlaneXY : public Shape
{
public:
	PlaneXY(float zCoord, glm::vec3 n, int materialIdx, bool flipNormal = false)
		: Shape(materialIdx), z(zCoord), normal(n)
	{
		if (flipNormal)
			normal *= -1;
	}
	bool Intersect(Ray& ray, float& tHit, SurfaceInteraction& intersection) const;


private:
	glm::vec3 normal;
	float z;
};

class PlaneXZ : public Shape
{
public:
	PlaneXZ(float yCoord, glm::vec3 n, int materialIdx, bool flipNormal = false)
		: Shape(materialIdx), y(yCoord), normal(n)
	{
		if (flipNormal)
			normal *= -1;
	}
	bool Intersect(Ray& ray, float& tHit, SurfaceInteraction& intersection) const;


private:
	glm::vec3 normal;
	float y;
};

class PlaneYZ : public Shape
{
public:
	PlaneYZ(float xCoord, glm::vec3 n, int materialIdx, bool flipNormal = false)
		: Shape(materialIdx), x(xCoord), normal(n)
	{
		if (flipNormal)
			normal *= -1;
	}
	bool Intersect(Ray& ray, float& tHit, SurfaceInteraction& intersection) const;


private:
	glm::vec3 normal;
	float x;
};