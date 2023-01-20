#include "PathTracer.hpp"
#include "../Utils/Utils.hpp"

#include <iostream>
#include <execution>

#include <glm/gtx/string_cast.hpp>
void PathTracer::Render(const Scene& scene, const Camera& camera)
{
	m_ActiveScene = &scene;
	if (frameIndex == 1)
		memset(accumulator, 0, m_FinalImage->GetWidth() * m_FinalImage->GetHeight() * sizeof(glm::vec3));

	
	//float invWidth = 1.f / (float)m_FinalImage->GetWidth();
	//float invHeight = 1.f / (float)m_FinalImage->GetHeight();

	m_ActiveCamera = &camera;
	std::for_each(std::execution::par, m_VerticalIter.begin(), m_VerticalIter.end(),
		[this](uint32_t y) {
			for (int x = 0; x < m_FinalImage->GetWidth(); x++)
			{
				//m_ActiveCamera = &camera;
				Ray ray;
				glm::vec2 randomOffset(Walnut::Random::Float(), Walnut::Random::Float());
				ray = m_ActiveCamera->getPrimaryRay(x, y, randomOffset);

				glm::vec3 color(TraceRay(ray,0));
				accumulator[x + y * m_FinalImage->GetWidth()] += color;

				glm::vec4 accumulated(glm::sqrt(accumulator[x + y * m_FinalImage->GetWidth()] / (float)frameIndex), 1.f);
				accumulated = glm::clamp(accumulated, glm::vec4(0.f), glm::vec4(1.f));

				m_ImageData[x + y * m_FinalImage->GetWidth()] = Utils::ConvertToRGBA(accumulated);
			}
		});
	//#pragma omp parallel for schedule(dynamic)
	//for (int y = 0; y < m_FinalImage->GetHeight(); y++)
	//{
	//	for (int x = 0; x < m_FinalImage->GetWidth(); x++)
	//	{
	//		m_ActiveCamera = &camera;
	//		Ray ray;
	//		glm::vec2 randomOffset(Walnut::Random::Float(), Walnut::Random::Float());
	//		ray = m_ActiveCamera->getPrimaryRay(x, y, randomOffset);

	//		glm::vec3 color(TraceRay(ray,0));
	//		accumulator[x + y * m_FinalImage->GetWidth()] += color;

	//		glm::vec4 accumulated(glm::pow(accumulator[x + y * m_FinalImage->GetWidth()] / (float)frameIndex, glm::vec3(gamma)), 1.f);
	//		accumulated = glm::clamp(accumulated, glm::vec4(0.f), glm::vec4(1.f));

	//		m_ImageData[x + y * m_FinalImage->GetWidth()] = Utils::ConvertToRGBA(accumulated);
	//	}
	//}
	m_FinalImage->SetData(m_ImageData);
	if (settings.Accumulate)
		frameIndex++;
	else
		Reset();
}
glm::vec3 PathTracer::TraceRay(Ray& ray, int depth)
{
	size_t maxBounces = 300;
	glm::vec3 throughput(1.f);
	for (size_t i = 0; i < maxBounces; i++)
	{
		Intersection isect;
		if (m_ActiveScene->Intersect(ray, isect))
		{
			SurfaceInteraction interaction = m_ActiveScene->getSurfaceProperties(ray, isect);
			Frame surf(interaction.hit_normal);

			Material mat = m_ActiveScene->materials[interaction.materialIdx];
			glm::vec3 hitColor = mat.albedo;
			float rrProb = 1;
			if (depth > 3)
			{
				rrProb = std::max({ hitColor.x, hitColor.y, hitColor.z });
				float rr = Walnut::Random::Float();
				if (rr > rrProb)
				{
					return glm::vec3(0.f);
				}
			}
			if (mat.textureIndex >= 0)
			{
				hitColor *= m_ActiveScene->textures[mat.textureIndex]->sampleImageTexture(interaction.uv.x, interaction.uv.y);
			}
			if (mat.radiance > 0.f)
			{
				return throughput * hitColor * mat.radiance;
			}
			if (mat.mirror)
			{
				glm::vec3 dir = surf.ToLocal(-ray.direction);
				dir.x *= -1;
				dir.y *= -1;
				dir = surf.ToWorld(dir);
				ray = getRay(rayPnt(ray, isect.t_hit) + interaction.hit_normal * 0.0001f, dir);
				throughput *= hitColor;
				continue;
			}
			else if (mat.glass)
			{
				float fres = fresnel(ray.direction, interaction.hit_normal, mat.ior);
				if (Walnut::Random::Float() < fres)
				{
					glm::vec3 dir = surf.ToLocal(-ray.direction);
					dir.x *= -1;
					dir.y *= -1;
					dir = surf.ToWorld(dir);
					if (glm::dot(interaction.hit_normal, dir) < 0)
					{
						ray = getRay(rayPnt(ray, isect.t_hit) - interaction.hit_normal * 0.0001f, glm::normalize(dir));
					}
					else
					{
						ray = getRay(rayPnt(ray, isect.t_hit) + interaction.hit_normal * 0.0001f, glm::normalize(dir));
					}
					throughput *= hitColor;
					continue;
				}
				else
				{
					glm::vec3 dir = refract(ray.direction, interaction.hit_normal, mat.ior);
					if (glm::dot(interaction.hit_normal, dir) < 0)
					{
						ray = getRay(rayPnt(ray, isect.t_hit) - interaction.hit_normal * 0.0001f, glm::normalize(dir));
					}
					else
					{
						ray = getRay(rayPnt(ray, isect.t_hit) + interaction.hit_normal * 0.0001f, glm::normalize(dir));
					}
					throughput *= hitColor;
					continue;
				}
			}
			else
			{
				float pdf;
				glm::vec3 dir = sampleCosineHemisphere(pdf);
				dir = surf.ToWorld(dir);
				glm::vec3 diffuse = glm::one_over_pi<float>() * hitColor * glm::dot(interaction.hit_normal, dir) / (pdf * rrProb);
				ray = getRay(rayPnt(ray, isect.t_hit) + interaction.hit_normal * 0.0001f, glm::normalize(dir));
				throughput *= diffuse;
				continue;
			}
		}
		else
		{
			return throughput * m_ActiveScene->getSkyColor(ray);
		}
	}
	return glm::vec3(0.f);
}

void PathTracer::OnResize(uint32_t width, uint32_t height)
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

	m_HorizontalIter.resize(width);
	m_VerticalIter.resize(height);
	for (uint32_t i = 0; i < width; i++)
		m_HorizontalIter[i] = i;
	for (uint32_t i = 0; i < height; i++)
		m_VerticalIter[i] = i;
}

