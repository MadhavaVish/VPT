#include "PathTracer.hpp"
#include "../Utils/Utils.hpp"

#include <iostream>

#include <glm/gtx/string_cast.hpp>
void PathTracer::Render(const Scene& scene, const Camera& camera)
{
	m_ActiveScene = &scene;
	if (frameIndex == 1)
		memset(accumulator, 0, m_FinalImage->GetWidth() * m_FinalImage->GetHeight() * sizeof(glm::vec3));
	float gamma = 1/2.2;
	float invWidth = 1.f / (float)m_FinalImage->GetWidth();
	float invHeight = 1.f / (float)m_FinalImage->GetHeight();
	float lensRadius = 0.008f;
	#pragma omp parallel for schedule(dynamic)
	for (int y = 0; y < m_FinalImage->GetHeight(); y++)
	{
		#pragma omp parallel for schedule(dynamic)
		for (int x = 0; x < m_FinalImage->GetWidth(); x++)
		{
			m_ActiveCamera = &camera;
			//Frame cam(m_ActiveCamera->GetDirection());
			Ray ray;
			//glm::vec2 rad = SampleConcentricDisc()*lensRadius;
			//glm::vec3 lol = cam.ToLocal(glm::vec3(rad, 0.f));
			ray.origin = m_ActiveCamera->GetPosition();
			//ray.origin = cam.ToWorld(cam.ToLocal(ray.origin) + lol);
			glm::vec2 randomOffset(Walnut::Random::Float(), Walnut::Random::Float());
			ray.direction = m_ActiveCamera->getPrimaryRay(x, y, randomOffset);
			//ray.direction = cam.ToWorld(cam.ToLocal(ray.direction) - lol);
			glm::vec3 color(TraceRay(ray,0));
			accumulator[x + y * m_FinalImage->GetWidth()] += color;
			glm::vec4 accumulated(glm::pow(accumulator[x + y * m_FinalImage->GetWidth()] / (float)frameIndex, glm::vec3(gamma)), 1.f);
			accumulated = glm::clamp(accumulated, glm::vec4(0.f), glm::vec4(1.f)) ;
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
static int maxBounces = 30;
glm::vec3 PathTracer::TraceRay(Ray& ray, int depth)
{
	
	if (depth == maxBounces)
	{
		return glm::vec3(0.f);
	}
	float closest = std::numeric_limits<float>::infinity();
	glm::vec3 color(0.f);
	SurfaceInteraction intersection;
	if (m_ActiveScene->Intersect(ray, closest, intersection))
	{
		
		Frame isect(intersection.hit_normal);
		Material mat = m_ActiveScene->materials[intersection.material];
		glm::vec3 hitColor = mat.albedo;
		float rrProb = 1;
		if (depth > 3)
		{
			rrProb = std::max({ hitColor.x, hitColor.y, hitColor.z });
			//rrProb = 0.5;
			float rr = Walnut::Random::Float();
			//rrProb = glm::min(rrProb, 0.95f);
			if (rr > rrProb)
			{
				return glm::vec3(0.f);
			}
		}
		if (mat.textureIndex >= 0)
		{
			hitColor *=m_ActiveScene->textures[mat.textureIndex]->sampleImageTexture(intersection.uv.x, intersection.uv.y);
		}
		if (mat.radiance > 0.f)
		{
			return hitColor *mat.radiance;
		}
		if (mat.mirror)
		{
			glm::vec3 dir = isect.ToLocal(-ray.direction);
			dir.x *= -1;
			dir.y *= -1;
			dir = isect.ToWorld(dir);
			Ray bounce;
			bounce.origin = ray(closest) + intersection.hit_normal*0.0001f;
			bounce.direction = glm::normalize(dir);
			return hitColor * TraceRay(bounce, depth+1);
		}
		else if (mat.glass)
		{
			float fres = fresnel(ray.direction, intersection.hit_normal, mat.ior);
			if (Walnut::Random::Float() < fres)
			{
				glm::vec3 dir = isect.ToLocal(-ray.direction);
				dir.x *= -1;
				dir.y *= -1;
				dir = isect.ToWorld(dir);
				Ray bounce;
				if (glm::dot(intersection.hit_normal, dir) < 0)
				{
					bounce.origin = ray(closest) - intersection.hit_normal * 0.0001f;
				}
				else
				{
					bounce.origin = ray(closest) + intersection.hit_normal * 0.0001f;
				}
				bounce.direction = glm::normalize(dir);
				return hitColor * TraceRay(bounce, depth + 1);
			}
			else
			{
				glm::vec3 dir = refract(ray.direction, intersection.hit_normal, mat.ior);
				Ray bounce;
				if (glm::dot(intersection.hit_normal, dir) < 0)
				{
					bounce.origin = ray(closest) - intersection.hit_normal * 0.0001f;
				}
				else
				{
					bounce.origin = ray(closest) + intersection.hit_normal * 0.0001f;
				}
				bounce.direction = glm::normalize(dir);
				return hitColor * TraceRay(bounce, depth + 1);
			}
		}
		else
		{
			//float pdf;
			//glm::vec3 dir = SampleUniformHemisphere(pdf);
			//dir = isect.ToWorld(dir);
			//Ray bounce;
			//bounce.origin = ray(closest) + intersection.hit_normal * 0.0001f;
			//bounce.direction = glm::normalize(dir);
			//return 2.f * hitColor * glm::dot(intersection.hit_normal, dir) * TraceRay(bounce, depth + 1);
			float pdf;
			glm::vec3 dir = sampleCosineHemisphere(pdf);
			dir = isect.ToWorld(dir);
			Ray bounce;
			bounce.origin = ray(closest) + intersection.hit_normal * 0.0001f;
			bounce.direction = glm::normalize(dir);
			glm:: vec3 throughput = glm::one_over_pi<float>() * hitColor * glm::dot(intersection.hit_normal, dir) / (pdf * rrProb);
			return TraceRay(bounce, depth + 1) * throughput;
		}
	}
	else
	{
		return m_ActiveScene->getSkyColor(ray);
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

