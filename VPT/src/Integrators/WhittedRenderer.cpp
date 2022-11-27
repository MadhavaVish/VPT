#include "WhittedRenderer.hpp"

#include <glm/gtc/constants.hpp>
#include <Walnut/Random.h>
#include <iostream>

namespace Utils {
	static uint32_t ConvertToRGBA(const glm::vec4& color) 
	{
		uint8_t r = (uint8_t)(color.r * 255.0f);
		uint8_t g = (uint8_t)(color.g * 255.0f);
		uint8_t b = (uint8_t)(color.b * 255.0f);
		uint8_t a = (uint8_t)(color.a * 255.0f);

		uint32_t result = (a << 24) | (b << 16) | (g << 8) | r;
		return result;
	}
}
void Whitted::Render(const Scene &scene, const Camera& camera)
{
	m_ActiveCamera = &camera;
	m_ActiveScene = &scene;
	#pragma omp parallel for
	for (int y = 0; y < m_FinalImage->GetHeight(); y++)
	{
		for (uint32_t x = 0; x < m_FinalImage->GetWidth(); x++)
		{
			glm::vec4 color = TraceRay(x,y);

			m_ImageData[x +y*m_FinalImage->GetWidth()] = Utils::ConvertToRGBA(color);

		}
	}
	m_FinalImage->SetData(m_ImageData);
}

glm::vec4 Whitted::TraceRay(int x, int y)
{
	Ray ray;
	ray.origin = m_ActiveCamera->GetPosition();
	ray.direction = m_ActiveCamera->GetRayDirections()[x + y * m_FinalImage->GetWidth()];
	glm::vec3 pointLight(1.f, 1.f, 1.f);
	float intensity = 3.f;
	float closest = std::numeric_limits<float>::infinity();
	int objIndex = m_ActiveScene->Intersect(ray, &closest);
	if (objIndex < 0)
	{
		return glm::vec4(0.1, 0.9, 1, 1);
	}
	else
	{
		glm::vec3 albedo = m_ActiveScene->spheres[objIndex].GetAlbedo(ray(closest));
		glm::vec3 normal = m_ActiveScene->spheres[objIndex].GetNormal(ray(closest));
		glm::vec3 directionToLight = glm::normalize(pointLight - ray(closest));
		float cosine = glm::dot(directionToLight, normal)*glm::one_over_pi<float>();
		glm::vec4 color(intensity * albedo * cosine, 1.f);
		color = glm::clamp(color, 0.f, 1.f);

		return color;
	}

}
void Whitted::OnResize(uint32_t width, uint32_t height)
{
	if (m_FinalImage)
	{
		if (m_FinalImage->GetWidth() == width && m_FinalImage->GetHeight() == height)
		{
			return;
		}
		m_FinalImage->Resize(width, height);
	}
	else
	{
		m_FinalImage = std::make_shared<Walnut::Image>(width, height, Walnut::ImageFormat::RGBA);
	}
	delete[] m_ImageData;
	m_ImageData = new uint32_t[width * height];


}