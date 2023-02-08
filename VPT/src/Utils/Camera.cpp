#include "Camera.hpp"
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>

#include <Walnut/Input/Input.h>
#include "../Scene/Materials/Material.hpp"
using namespace Walnut;

Camera::Camera(float fov) :m_VerticalFOV(fov)
{
	forwardDir = glm::vec3(0.f, 0.f, -1.f);
	position = glm::vec3(0.f, 0.f, 8.f);
}
bool Camera::OnUpdate(float ts)
{
	glm::vec2 mousePos = Input::GetMousePosition();
	glm::vec2 delta = (mousePos - m_LastMousePosition) * 0.002f;
	m_LastMousePosition = mousePos;

	if (!Input::IsMouseButtonDown(MouseButton::Right))
	{
		Input::SetCursorMode(CursorMode::Normal);
		return false;
	}

	Input::SetCursorMode(CursorMode::Locked);

	bool moved = false;

	constexpr glm::vec3 upDirection(0.0f, 1.0f, 0.0f);
	glm::vec3 rightDirection = glm::cross(forwardDir, upDirection);

	float speed = 5.0f;

	// Movement
	if (Input::IsKeyDown(KeyCode::W))
	{
		position += forwardDir * speed * ts;
		moved = true;
	}
	else if (Input::IsKeyDown(KeyCode::S))
	{
		position -= forwardDir * speed * ts;
		moved = true;
	}
	if (Input::IsKeyDown(KeyCode::A))
	{
		position -= rightDirection * speed * ts;
		moved = true;
	}
	else if (Input::IsKeyDown(KeyCode::D))
	{
		position += rightDirection * speed * ts;
		moved = true;
	}
	if (Input::IsKeyDown(KeyCode::Q))
	{
		position -= _y * speed * ts;
		moved = true;
	}
	else if (Input::IsKeyDown(KeyCode::E))
	{
		position += _y * speed * ts;
		moved = true;
	}

	// Rotation
	if (delta.x != 0.0f || delta.y != 0.0f)
	{
		float pitchDelta = delta.y * GetRotationSpeed();
		float yawDelta = delta.x * -GetRotationSpeed();

		glm::quat q = glm::normalize(glm::cross(glm::angleAxis(-pitchDelta, rightDirection),
			glm::angleAxis(-yawDelta, _y)));
		forwardDir = glm::rotate(q, forwardDir);

		moved = true;
	}

	if (moved)
	{
		RecalculateView();
	}
	return moved;
}

void Camera::OnResize(uint32_t width, uint32_t height)
{
	if (width == m_ViewportWidth && height == m_ViewportHeight)
	{
		RecalculateView();
		return;
	}
	m_ViewportWidth = width;
	m_ViewportHeight = height;
	invHeight = 1.f / (float)m_ViewportHeight;
	invWidth = 1.f / (float)m_ViewportWidth;
	aspect = (float)m_ViewportHeight/ (float)m_ViewportWidth;
	RecalculateView();
}
const Ray Camera::getPrimaryRay(const int x, const int y, const glm::vec2 offset) const
{
	float u = (x + offset.x) * invWidth;
	float v = (y + offset.y) * invHeight;
	u = u * 2 - 1;
	v = v * 2 - 1;
	glm::vec3 defocus{ 0.f };
	if (settings.f_stop < 30.f)
	{
		//glm::vec2 rad{ Walnut::Random::Float()-0.5f, Walnut::Random::Float() - 0.5f };
		//rad *= aperture;
		glm::vec2 rad = SampleConcentricDisc() * aperture;
		defocus = _x * rad.x + _y * rad.y;

	}
	glm::vec3 dir = glm::normalize(glm::normalize(forwardDir * settings.focus_dist + u * horizontal - v * vertical) * settings.focus_dist - defocus);

	return getRay(position+defocus, dir);
}
void Camera::RecalculateView()
{
	float w = settings.sensor_width*0.001f / (2 * settings.focal_length * 0.001f);
	aperture = 0.5f * settings.focal_length * 0.001f / settings.f_stop;
	float width = 2 * w;
	float height = aspect * width;

	_z = forwardDir;
	_x = glm::normalize(glm::cross(_z, { 0.f, 1.f, 0.f }));
	_y = glm::cross(_z, _x);

	horizontal = settings.focus_dist * width * _x/2.f;
	vertical = settings.focus_dist * height * _y/2.f;
}