#include "PathTracer.hpp"
#include "../Utils/Utils.hpp"

#include <iostream>

#include <glm/gtx/string_cast.hpp>
void PathTracer::Render(const Scene& scene, const Camera& camera)
{
	m_ActiveScene = &scene;
	if (frameIndex == 1)
		memset(accumulator, 0, m_FinalImage->GetWidth() * m_FinalImage->GetHeight() * sizeof(glm::vec3));

	float gamma = 1/2.2f;
	float invWidth = 1.f / (float)m_FinalImage->GetWidth();
	float invHeight = 1.f / (float)m_FinalImage->GetHeight();
	float lensRadius = 0.008f;
	#pragma omp parallel for schedule(dynamic)
	for (int y = 0; y < m_FinalImage->GetHeight(); y++)
	{
		for (int x = 0; x < m_FinalImage->GetWidth(); x++)
		{
			m_ActiveCamera = &camera;
			Ray ray;
			glm::vec2 randomOffset(Walnut::Random::Float(), Walnut::Random::Float());
			ray = m_ActiveCamera->getPrimaryRay(x, y, randomOffset);

			glm::vec3 color(TraceRay(ray,0));
			accumulator[x + y * m_FinalImage->GetWidth()] += color;

			glm::vec4 accumulated(glm::pow(accumulator[x + y * m_FinalImage->GetWidth()] / (float)frameIndex, glm::vec3(gamma)), 1.f);
			accumulated = glm::clamp(accumulated, glm::vec4(0.f), glm::vec4(1.f));

			if (settings.Vignette)
			{
				float exp = 1 / settings.VignetteAmount;
				float _x = x * invWidth, _y = y * invHeight;
				accumulated *= (1 + glm::pow(8 * _x * _y * (1 - _x) * (1 - _y), exp)) / 2.f;

			}
			m_ImageData[x + y * m_FinalImage->GetWidth()] = Utils::ConvertToRGBA(accumulated);
		}
	}
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
		if (depth == maxBounces)
		{
			return glm::vec3(0.f);
		}
		//float closest = std::numeric_limits<float>::infinity();
		glm::vec3 color(0.f);
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
				Ray bounce;
				bounce.origin = ray(isect.t_hit) + interaction.hit_normal * 0.0001f;
				bounce.direction = glm::normalize(dir);
				ray = bounce;
				//return hitColor * TraceRay(bounce, depth + 1);
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
					Ray bounce;
					if (glm::dot(interaction.hit_normal, dir) < 0)
					{
						bounce.origin = ray(isect.t_hit) - interaction.hit_normal * 0.0001f;
					}
					else
					{
						bounce.origin = ray(isect.t_hit) + interaction.hit_normal * 0.0001f;
					}
					bounce.direction = glm::normalize(dir);
					ray = bounce;
					//return hitColor * TraceRay(bounce, depth + 1);
					throughput *= hitColor;
					continue;
				}
				else
				{
					glm::vec3 dir = refract(ray.direction, interaction.hit_normal, mat.ior);
					Ray bounce;
					if (glm::dot(interaction.hit_normal, dir) < 0)
					{
						bounce.origin = ray(isect.t_hit) - interaction.hit_normal * 0.0001f;
					}
					else
					{
						bounce.origin = ray(isect.t_hit) + interaction.hit_normal * 0.0001f;
					}
					bounce.direction = glm::normalize(dir);
					ray = bounce;
					//return hitColor * TraceRay(bounce, depth + 1);
					throughput *= hitColor;
					continue;
				}
			}
			else
			{
				float pdf;
				glm::vec3 dir = sampleCosineHemisphere(pdf);
				dir = surf.ToWorld(dir);
				Ray bounce;
				bounce.origin = ray(isect.t_hit) + interaction.hit_normal * 0.0001f;
				bounce.direction = glm::normalize(dir);
				glm::vec3 diffuse = glm::one_over_pi<float>() * hitColor * glm::dot(interaction.hit_normal, dir) / (pdf * rrProb);
				//return TraceRay(bounce, depth + 1) * throughput;
				ray = bounce;
				throughput *= diffuse;
				continue;
			}
		}
		else
		{
			return throughput * m_ActiveScene->getSkyColor(ray);
		}
	}
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
}

