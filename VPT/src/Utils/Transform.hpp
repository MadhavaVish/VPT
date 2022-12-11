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

struct Frame
{
	Frame(const glm::vec3 up)
	{
		glm::vec3 z = m_Z = glm::normalize(up);
		glm::vec3 temp = glm::abs(z.x) > 0.99f ? glm::vec3(0.f, 1.f, 0.f) : glm::vec3(1.f, 0.f, 0.f);
		m_Y = glm::normalize(glm::cross(m_Z, temp));
		m_X = glm::cross(m_Y, m_Z);
	}
	void setFromUp(glm::vec3 up)
	{
		glm::vec3 z = m_Z = glm::normalize(up);
		glm::vec3 temp = glm::abs(z.x) > 0.99f ? glm::vec3(0.f, 1.f, 0.f) : glm::vec3(1.f, 0.f, 0.f);
		m_Y = glm::normalize(glm::cross(m_Z, temp));
		m_X = glm::cross(m_Y, m_Z);
	}
	glm::vec3 ToWorld(const glm::vec3& a) const
	{
		return m_X * a.x + m_Y * a.y + m_Z * a.z;
	}
	glm::vec3 ToLocal(const glm::vec3 a) const
	{
		return glm::vec3(glm::dot(a, m_X), glm::dot(a, m_Y), glm::dot(a, m_Z));
	}
	glm::vec3 m_X, m_Y, m_Z;
};