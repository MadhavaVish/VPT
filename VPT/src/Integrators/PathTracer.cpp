#include "PathTracer.hpp"
#include "../Utils/Utils.hpp"

#include <iostream>

#include <glm/gtx/string_cast.hpp>
void PathTracer::Render(const Scene& scene, const Camera& camera)
{
	m_ActiveScene = &scene;
	if (frameIndex == 1)
		memset(accumulator, 0, m_FinalImage->GetWidth() * m_FinalImage->GetHeight() * sizeof(glm::vec3));
	float gamma = 1 / 2.2f;
	#pragma omp parallel for schedule(dynamic)
	for (int y = 0; y < m_FinalImage->GetHeight(); y++)
	{
		m_ActiveCamera = &camera;
		for (uint32_t x = 0; x < m_FinalImage->GetWidth(); x++)
		{
			Ray ray;
			ray.origin = m_ActiveCamera->GetPosition();
			ray.direction = m_ActiveCamera->GetRayDirections()[x + y * m_FinalImage->GetWidth()];
			glm::vec3 color(TraceRay(ray));
			accumulator[x + y * m_FinalImage->GetWidth()] += color;

			glm::vec4 accumulated(glm::pow(accumulator[x + y * m_FinalImage->GetWidth()] / (float)frameIndex, glm::vec3(gamma)), 1.f);
			//accumulated = glm::pow(accumulated, glm::vec3(gamma));
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

static void sampleBRDF(const Material& mat, glm::vec3& incoming, glm::vec3& outgoing, float& pdf)
{
	
		outgoing = sampleCosineHemisphere(pdf);
}
glm::vec3 evalBRDF(const Material& mat, glm::vec3& incoming, glm::vec3& outgoing)
{
	glm::vec3 diffuseComponent = mat.albedo * glm::one_over_pi<float>();

	return diffuseComponent;
}
glm::vec3 PathTracer::TraceRay(Ray& ray)
{
	int maxBounces = 7;
	float closest = std::numeric_limits<float>::infinity();
	glm::vec3 color(0.f);
	SurfaceInteraction intersection;
	if (m_ActiveScene->Intersect(ray, closest, intersection))
	{
		glm::vec3 outgoing;
		float prob;
		Frame hitFrame(intersection.hit_normal);
		sampleBRDF(m_ActiveScene->materials[intersection.material], hitFrame.ToLocal(ray.direction), outgoing, prob);
		if (glm::length(outgoing) < 0.0001)
		{
			return glm::vec3(0.f);
		}

		color = glm::vec3(0.f);
		//color += m_ActiveScene->materials[intersection.material].albedo * glm::one_over_pi<float>() * glm::dot(intersection.hit_normal, outgoing) / prob;
		color += evalBRDF(m_ActiveScene->materials[intersection.material], hitFrame.ToLocal(ray.direction), outgoing) / prob;
		outgoing = hitFrame.ToWorld(outgoing);
		color *= glm::dot(intersection.hit_normal, outgoing);
		ray.origin = ray(closest) + intersection.hit_normal * 0.0001f;
		ray.direction = outgoing;
		for (int bounce = 0; bounce < maxBounces; bounce++)
		{
			if (bounce == (maxBounces -1))
			{
				color *= glm::vec3(0.f);
				break;
			}
			if (m_ActiveScene->Intersect(ray, closest, intersection))
			{
				glm::vec3 outgoing;
				float prob;
				hitFrame.setFromUp(intersection.hit_normal);
				sampleBRDF(m_ActiveScene->materials[intersection.material], hitFrame.ToLocal(ray.direction), outgoing, prob);
				if (glm::length(outgoing) < 0.0001)
				{
					return glm::vec3(0.f);
				}
				//glm::vec3 f = glm::vec3(0.f);
				//f += m_ActiveScene->materials[intersection.material].albedo * glm::one_over_pi<float>() * glm::dot(intersection.hit_normal, outgoing) / prob;
				glm::vec3 f = evalBRDF(m_ActiveScene->materials[intersection.material], hitFrame.ToLocal(ray.direction), outgoing) / prob;

				outgoing = hitFrame.ToWorld(outgoing);
				f *= glm::dot(intersection.hit_normal, outgoing);
				color *= f;
				ray.origin = ray(closest) + intersection.hit_normal * 0.0001f;
				ray.direction = outgoing;
			}
			else
			{
				color*= glm::vec3(1.f);
				break;
			}
		}
	}
	else
	{
		return glm::vec3(1.f);
	}
	return color;
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

