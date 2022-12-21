#include "Visualizer.hpp"
#include "../Utils/Utils.hpp"

#include <iostream>

#include <glm/gtx/string_cast.hpp>
void Visualizer::Render(const Scene& scene, const Camera& camera)
{
	m_ActiveScene = &scene;
	if (frameIndex == 1)
		memset(accumulator, 0, m_FinalImage->GetWidth() * m_FinalImage->GetHeight() * sizeof(glm::vec3));

	float gamma = 1 / 2.2f;
	#pragma omp parallel for schedule(dynamic)
	for (int y = 0; y < m_FinalImage->GetHeight(); y++)
	{
		for (uint32_t x = 0; x < m_FinalImage->GetWidth(); x++)
		{
			m_ActiveCamera = &camera;
			Ray ray;
			glm::vec2 randomOffset(Walnut::Random::Float(), Walnut::Random::Float());
			ray = m_ActiveCamera->getPrimaryRay(x, y, randomOffset);

			glm::vec3 color(TraceRay(ray, 0));
			accumulator[x + y * m_FinalImage->GetWidth()] += color;

			glm::vec4 accumulated(glm::pow(accumulator[x + y * m_FinalImage->GetWidth()] / (float)frameIndex, glm::vec3(gamma)), 1.f);
			accumulated = glm::clamp(accumulated, glm::vec4(0.f), glm::vec4(1.f));

			m_ImageData[x + y * m_FinalImage->GetWidth()] = Utils::ConvertToRGBA(accumulated);
		}
	}
	m_FinalImage->SetData(m_ImageData);
	if (settings.Accumulate)
		frameIndex++;
	else
		Reset();
}
static float near = 0.1f;
static float far = 3000.f;

glm::vec3 Visualizer::TraceRay(Ray& ray, int depth)
{
	Intersection isect;
	if (m_ActiveScene->Intersect(ray, isect))
	{
		SurfaceInteraction interaction = m_ActiveScene->getSurfaceProperties(ray, isect);
		Frame surf(interaction.hit_normal);

		Material mat = m_ActiveScene->materials[interaction.materialIdx];
		glm::vec3 hitColor = mat.albedo;
		if (mat.textureIndex >= 0)
		{
			hitColor *= m_ActiveScene->textures[mat.textureIndex]->sampleImageTexture(interaction.uv.x, interaction.uv.y);
		}
		return hitColor;
		/*glm::vec3 normalColor = interaction.hit_normal * 0.5f + 0.5f;
		return normalColor;*/
		//float d = far / (far - near) + (1.f / isect.t_hit) * (-far * near) / (far - near);
		/*glm::vec3 depth(isect.t_hit/10.f);
		return depth;*/
	}
	else
	{
		//return glm::vec3(0.f);
		return m_ActiveScene->getSkyColor(ray);
	}
	
	return glm::vec3(0.f);
}

void Visualizer::OnResize(uint32_t width, uint32_t height)
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
	delete[] accumulator;
	accumulator = new glm::vec3[width * height];
}

