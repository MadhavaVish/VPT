#pragma once
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>
#include <iostream>
#include "Ray.hpp"
class Camera
{
public:
	struct Settings
	{
		float f_stop{ 32.f };
		float focus_dist{ 1.f };
		float focal_length{ 28.f };
		float sensor_width{ 36.f };
	};

	Camera() : Camera(16.1f){};
	Camera(float fov);
		
	bool OnUpdate(float ts);
	void OnResize(uint32_t width, uint32_t height);
	const Ray getPrimaryRay(const int x, const int y, const glm::vec2 offset = glm::vec2(0.f)) const;
	float GetRotationSpeed()
	{
		return 0.3f;
	}

	const glm::vec3& GetPosition() const { return position; }
	const glm::vec3& GetDirection() const { return forwardDir; }
	Settings settings;
private:
	void RecalculateView();
private:
	glm::vec3 position{ 0.0f, 0.0f, 0.0f };
	glm::vec3 forwardDir{ 0.0f, 0.0f, 0.0f };
	glm::vec3 horizontal{ 0.0f, 0.0f, 0.0f }, vertical{ 0.0f, 0.0f, 0.0f };
	glm::vec3 _x{ 0.0f, 0.0f, 0.0f }, _y{ 0.0f, 0.0f, 0.0f }, _z{ 0.0f, 0.0f, 0.0f };
	float m_VerticalFOV = 45.0f;
	float aspect{ 0.f };
	float aperture{ 1.f };
	glm::vec2 m_LastMousePosition{ 0.0f, 0.0f };

	uint32_t m_ViewportWidth = 0, m_ViewportHeight = 0;

};