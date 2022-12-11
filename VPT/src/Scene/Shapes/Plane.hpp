#pragma once
#include "Shape.hpp"


class Plane : public Shape
{
public:

	Plane(glm::vec3 n, float d, int materialIdx) ;
	bool Intersect(Ray& ray, float& tHit, SurfaceInteraction &intersection) const;
private:
	Frame local;
	glm::vec3 normal;
	float distance;
};

class PlaneXY : public Shape
{
public:		
	PlaneXY(float zCoord, glm::vec3 n, int materialIdx, bool flipNormal = false)
		: Shape(materialIdx), z(-zCoord), normal(n), local(Frame(n))
	{
		if (flipNormal)
			normal *= -1;
	}
	bool Intersect(Ray& ray, float& tHit, SurfaceInteraction& intersection) const;


private:
	Frame local;
	glm::vec3 normal;
	float z;
};

class PlaneXZ : public Shape
{
public:
	PlaneXZ(float yCoord, glm::vec3 n, int materialIdx, bool flipNormal = false)
		: Shape(materialIdx), y(-yCoord), normal(n), local((n))
	{
		if (flipNormal)
			normal *= -1;
	}
	bool Intersect(Ray& ray, float& tHit, SurfaceInteraction& intersection) const;


private:
	Frame local;
	glm::vec3 normal;
	float y;
};

class PlaneYZ : public Shape
{
public:
	PlaneYZ(float xCoord, glm::vec3 n, int materialIdx, bool flipNormal = false)
		: Shape(materialIdx), x(-xCoord), normal(n), local(Frame(n))
	{
		if (flipNormal)
			normal *= -1;
	}
	bool Intersect(Ray& ray, float& tHit, SurfaceInteraction& intersection) const;


private:
	Frame local;
	glm::vec3 normal;
	float x;
};