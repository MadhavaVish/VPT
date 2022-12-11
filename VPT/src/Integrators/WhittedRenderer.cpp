#include "WhittedRenderer.hpp"
#include "../Utils/Utils.hpp"
#include <glm/gtc/constants.hpp>
#include <iostream>

void Whitted::Render(const Scene &scene, const Camera& camera)
{
	m_ActiveScene = &scene;
	float gamma = 1/2.2;
	const float invWidth = 1.f / m_FinalImage->GetWidth();
	const float invHeight = 1.f / m_FinalImage->GetWidth();
	#pragma omp parallel for schedule(dynamic)
	for (int y = 0; y < m_FinalImage->GetHeight(); y++)
	{
		m_ActiveCamera = &camera;
		for (uint32_t x = 0; x < m_FinalImage->GetWidth(); x++)
		{
			Ray ray;
			ray.origin = m_ActiveCamera->GetPosition();
			glm::vec2 randomOffset(Walnut::Random::Float() * invWidth, Walnut::Random::Float() * invHeight);
			ray.direction = m_ActiveCamera->getPrimaryRay(x, y, randomOffset);
			//ray.direction = m_ActiveCamera->GetRayDirections()[x + y * m_FinalImage->GetWidth()];
			//glm::vec4 color(TraceRay(ray, 0), 1.f);
			glm::vec4 color(glm::pow(TraceRay(ray, 0), glm::vec3(gamma)), 1.f);
			m_ImageData[x +y*m_FinalImage->GetWidth()] = Utils::ConvertToRGBA(color);
		}
	}
	m_FinalImage->SetData(m_ImageData);
}
constexpr int maxDepth = 50;

glm::vec3 Whitted::TraceRay(Ray &ray, int depth)
{
	if (depth > maxDepth)
	{
		return glm::vec3(0.4f);
	}
	glm::vec3 pointLight(1.f, 1.f, 1.f);
	float intensity = 3;
	
	float closest = std::numeric_limits<float>::infinity();
	SurfaceInteraction intersection;
	bool hit = m_ActiveScene->Intersect(ray, closest, intersection);

	if (!hit)
	{
		return m_ActiveScene->getSkyColor(ray);
	}
	else
	{
		Material mat = m_ActiveScene->materials[intersection.material];
		glm::vec3 hitColor = mat.albedo;
		if (mat.textureIndex >= 0)
		{
			hitColor *= m_ActiveScene->textures[mat.textureIndex]->sampleImageTexture(intersection.uv.x, intersection.uv.y);
		}
		glm::vec3 surfaceNormal = intersection.hit_normal;
		glm::vec3 hitPoint = ray(closest);
		if (m_ActiveScene->materials[intersection.material].glass)
		{
			glm::vec3 refl = glm::normalize(reflect(ray.direction, surfaceNormal));
			glm::vec3 refr = glm::normalize(refract(ray.direction, surfaceNormal, m_ActiveScene->materials[intersection.material].ior));
			glm::vec3 reflO = glm::dot(refl, surfaceNormal) < 0 ? hitPoint - surfaceNormal* 0.0001f : hitPoint + surfaceNormal* 0.0001f;
			glm::vec3 refrO = glm::dot(refr, surfaceNormal) < 0 ? hitPoint - surfaceNormal * 0.0001f : hitPoint + surfaceNormal * 0.0001f;
			Ray reflected(reflO, refl);
			Ray refracted(refrO, refr);
			float fr = fresnel(ray.direction, surfaceNormal, m_ActiveScene->materials[intersection.material].ior);
			if (glm::dot(ray.direction, surfaceNormal) > 0)
			{
				glm::vec3 path = hitPoint - ray.origin;
				float pathLength = glm::sqrt(glm::dot(path, path));
				glm::vec3 attentuation = glm::exp(-(hitColor) * 0.6f * pathLength);
				//glm::vec3 attentuation(1.f);
				glm::vec3 refrColor = hitColor * attentuation * TraceRay(refracted, depth + 1);
				glm::vec3 reflColor = hitColor * TraceRay(reflected, depth + 1);
				return (fr) * reflColor + (1-fr) * refrColor;
				//return refrColor;
			}
			else
			{
				glm::vec3 refrColor = hitColor*TraceRay(refracted, depth + 1);
				glm::vec3 reflColor = hitColor*TraceRay(reflected, depth + 1);
				return fr * reflColor + (1 - fr) * refrColor;
				//return refrColor;
			}

		}
		else if (m_ActiveScene->materials[intersection.material].mirror)
		{
			glm::vec3 refl = glm::normalize(reflect(ray.direction, surfaceNormal));
			glm::vec3 reflO = glm::dot(refl, surfaceNormal) < 0 ? hitPoint - surfaceNormal * 0.0001f : hitPoint + surfaceNormal * 0.0001f;
			Ray reflected(reflO, refl);
			return hitColor*TraceRay(reflected, depth + 1);
		}
		else
		{
			glm::vec3 directionToLight = pointLight-hitPoint;
			float dist2 = glm::dot(directionToLight, directionToLight);
			float dist = glm::sqrt(dist2);
			directionToLight = glm::normalize(directionToLight);
			Ray shadow(hitPoint+ surfaceNormal *0.0001f, directionToLight);
			float cosL = glm::dot(surfaceNormal, directionToLight);
			float occluded = std::numeric_limits<float>::infinity();
			SurfaceInteraction nop;
			if (m_ActiveScene->Intersect(shadow, occluded, nop))
			{
				//float dist = glm::sqrt(dist2);
				if (occluded*occluded < dist2) {
					return glm::vec3(0.f);
				}
			}
			glm::vec3 albedo = m_ActiveScene->materials[intersection.material].albedo;
			//float diffuse = m_ActiveScene->materials[intersection.material].roughness;
			//float diffuse = intensity * cosL * glm::one_over_pi<float>()/dist2;
			//float specular = 1 - diffuse;
			glm::vec3 diffColor = hitColor * intensity * cosL * glm::one_over_pi<float>()/ dist2;
			//glm::vec3 reflectionDir = reflect(-directionToLight, surfaceNormal);
			//if (specular < 0.00001f) {
				return glm::clamp(diffColor, 0.f, 1.f);
			//}
			//else
			//{
			//	//glm::vec3 reflectionDir = glm::normalize(reflect(ray.direction, surfaceNormal));
			//	float exponent = glm::max(2 / glm::pow(diffuse, 4.f)  - 2, 0.01f);
			//	float specColor = intensity * glm::pow(glm::max(0.f, glm::dot(reflectionDir, -ray.direction)), exponent);
			//	//glm::vec3 reflO = glm::dot(reflectionDir, surfaceNormal) < 0 ? hitPoint - surfaceNormal * 0.0001f : hitPoint + surfaceNormal * 0.0001f;
			//	//Ray reflected(reflO, reflectionDir);
			//	return glm::clamp(diffColor * diffuse + specular * specColor, 0.f, 1.f);
			//	
			//}
			
		}
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