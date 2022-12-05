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
	#pragma omp parallel for schedule(dynamic)
	for (int y = 0; y < m_FinalImage->GetHeight(); y++)
	{
		m_ActiveCamera = &camera;
		for (uint32_t x = 0; x < m_FinalImage->GetWidth(); x++)
		{
			Ray ray;
			ray.origin = m_ActiveCamera->GetPosition();
			ray.direction = m_ActiveCamera->GetRayDirections()[x + y * m_FinalImage->GetWidth()];
			glm::vec3 color(TraceRay(ray,0));
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
static float BaseReflectance(const float& eta1, const float& eta2)
{
	return glm::pow((eta2 - eta1) / (eta2 + eta1), 2);
}
//static float Fresnel(const float& baseReflectance, const float& cosTheta)
//{
//	return baseReflectance + (1 - baseReflectance) * glm::pow((1 - cosTheta), 5);
//}
static void sampleBRDF(const Material& mat, glm::vec3& incoming, glm::vec3& outgoing, float& pdf)
{
	if (mat.mirror)
	{
		glm::vec3 reflected(-incoming.x, -incoming.y, incoming.z);
		outgoing = glm::normalize(reflected);
		pdf = 1;
	}
	else if (mat.glass)
	{
		//float fresnel = Fresnel(BaseReflectance(mat.ior, 1), incoming.z);
		float fre = fresnel(incoming, glm::vec3(0.f, 0.f, 1.f), mat.ior);
		if (Walnut::Random::Float() < fre)
		{
			glm::vec3 reflected(-incoming.x, -incoming.y, incoming.z);
			outgoing = glm::normalize(reflected);
			pdf = 1;
		}
		else
		{
			outgoing = glm::normalize(refract(-incoming, glm::vec3(0.f, 0.f, 1.f), mat.ior));
			pdf = 1;
		}
	}
	else {
		outgoing = SampleUniformHemisphere(pdf);
	}
}
glm::vec3 evalBRDF(const Material& mat, glm::vec3& incoming, glm::vec3& outgoing)
{
	if (mat.mirror)
	{
		return mat.albedo;
	}
	else if (mat.glass)
	{
		float fre = fresnel(incoming, glm::vec3(0.f, 0.f, 1.f), mat.ior);
		return mat.albedo*(1-fre);
	}
	else
	{

		float kd = std::max({ mat.albedo.r, mat.albedo.g, mat.albedo.b });
		float ks = mat.ks;
		float norm = 1 / (kd + ks);
		kd *= norm;
		ks *= norm;

		glm::vec3 diffuseComponent = kd * mat.albedo * glm::one_over_pi<float>();
		glm::vec3 specular = ks * glm::vec3(1.f) * (mat.exponent + 2) * glm::one_over_two_pi<float>() * glm::pow(glm::dot(incoming, outgoing), mat.exponent);
		return diffuseComponent + specular;
	}
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
		if (mat.textureIndex >= 0)
		{
			hitColor *=m_ActiveScene->textures[mat.textureIndex]->sampleImageTexture(intersection.uv.x, intersection.uv.y);
		}
		if (mat.radiance > 0.f)
		{
			return mat.albedo*mat.radiance;
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
			/*float pdf;
			glm::vec3 dir = SampleUniformHemisphere(pdf);
			dir = isect.ToWorld(dir);
			Ray bounce;
			bounce.origin = ray(closest) + intersection.hit_normal * 0.0001f;
			bounce.direction = glm::normalize(dir);
			return 2.f * mat.albedo * glm::dot(intersection.hit_normal, dir) * TraceRay(bounce, depth + 1);*/
			float pdf;
			glm::vec3 dir = sampleCosineHemisphere(pdf);
			dir = isect.ToWorld(dir);
			Ray bounce;
			bounce.origin = ray(closest) + intersection.hit_normal * 0.0001f;
			bounce.direction = glm::normalize(dir);
			return glm::one_over_pi<float>()* hitColor * glm::dot(intersection.hit_normal, dir)* TraceRay(bounce, depth + 1) /pdf;
		}
	}
	else
	{
		//return glm::vec3(1.f);
		return m_ActiveScene->getSkyColor(ray);
	}
}
//glm::vec3 PathTracer::TraceRay(Ray& ray)
//{
//	int maxBounces = 7;
//	float closest = std::numeric_limits<float>::infinity();
//	glm::vec3 color(0.f);
//	SurfaceInteraction intersection;
//	if (m_ActiveScene->Intersect(ray, closest, intersection))
//	{
//		glm::vec3 outgoing;
//		float prob;
//		Frame hitFrame(intersection.hit_normal);
//		sampleBRDF(m_ActiveScene->materials[intersection.material], hitFrame.ToLocal(-ray.direction), outgoing, prob);
//		if (glm::length(outgoing) < 0.0001f)
//		{
//			return glm::vec3(0.f);
//		}
//
//		color = glm::vec3(0.f);
//		//color += m_ActiveScene->materials[intersection.material].albedo * glm::one_over_pi<float>() * glm::dot(intersection.hit_normal, outgoing) / prob;
//		color += evalBRDF(m_ActiveScene->materials[intersection.material], hitFrame.ToLocal(ray.direction), outgoing) / prob;
//		outgoing = hitFrame.ToWorld(outgoing);
//		glm::vec3 reflO = glm::dot(outgoing, intersection.hit_normal) < 0 ? ray(closest) - intersection.hit_normal * 0.0001f : ray(closest) + intersection.hit_normal * 0.0001f;
//		if (!m_ActiveScene->materials[intersection.material].glass && !m_ActiveScene->materials[intersection.material].mirror)
//		{
//			color *= glm::dot(intersection.hit_normal, outgoing);
//		}
//		//if(glm::dot(outgoing, intersection.hit_normal))
//		//ray.origin = ray(closest) + intersection.hit_normal * 0.0001f;
//		ray.origin = reflO;
//		ray.direction = outgoing;
//		for (int bounce = 0; bounce < maxBounces; bounce++)
//		{
//			if (bounce == (maxBounces -1))
//			{
//				color *= glm::vec3(0.f);
//				break;
//			}
//			if (m_ActiveScene->Intersect(ray, closest, intersection))
//			{
//				glm::vec3 outgoing;
//				float prob;
//				hitFrame.setFromUp(intersection.hit_normal);
//				sampleBRDF(m_ActiveScene->materials[intersection.material], hitFrame.ToLocal(-ray.direction), outgoing, prob);
//				/*if (glm::length(outgoing) < 0.0001f)
//				{
//					return glm::vec3(0.f);
//				}*/
//				//glm::vec3 f = glm::vec3(0.f);
//				//f += m_ActiveScene->materials[intersection.material].albedo * glm::one_over_pi<float>() * glm::dot(intersection.hit_normal, outgoing) / prob;
//				glm::vec3 f = evalBRDF(m_ActiveScene->materials[intersection.material], hitFrame.ToLocal(-ray.direction), outgoing) / prob;
//
//				outgoing = hitFrame.ToWorld(outgoing);
//				glm::vec3 reflO = glm::dot(outgoing, intersection.hit_normal) < 0 ? ray(closest) - intersection.hit_normal * 0.0001f : ray(closest) + intersection.hit_normal * 0.0001f;
//				if (!m_ActiveScene->materials[intersection.material].glass && !m_ActiveScene->materials[intersection.material].mirror)
//				{
//					color *= glm::dot(intersection.hit_normal, outgoing);
//				}
//				//f *= glm::dot(intersection.hit_normal, outgoing);
//				color *= f;
//				//ray.origin = ray(closest) + intersection.hit_normal * 0.0001f;
//				ray.origin = reflO;
//				ray.direction = outgoing;
//			}
//			else
//			{
//				color*= glm::vec3(1.f);
//				break;
//			}
//		}
//	}
//	else
//	{
//		return glm::vec3(1.f);
//	}
//	return color;
//}

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

