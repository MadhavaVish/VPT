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
	m_ActiveScene = &scene;
	#pragma omp parallel for schedule(dynamic)
	for (int y = 0; y < m_FinalImage->GetHeight(); y++)
	{
		m_ActiveCamera = &camera;
		for (uint32_t x = 0; x < m_FinalImage->GetWidth(); x++)
		{
			Ray ray;
			ray.origin = m_ActiveCamera->GetPosition();
			ray.direction = m_ActiveCamera->GetRayDirections()[x + y * m_FinalImage->GetWidth()];
			glm::vec4 color(TraceRay(ray, 0), 1.f);

			m_ImageData[x +y*m_FinalImage->GetWidth()] = Utils::ConvertToRGBA(color);
		}
	}
	m_FinalImage->SetData(m_ImageData);
}
constexpr int maxDepth = 5;
float fresnel(const glm::vec3& I, const glm::vec3& N, const float& ior)
{
	float cosi = glm::clamp(glm::dot (I, N), -1.f, 1.f);
	float etai = 1, etat = ior;
	if (cosi > 0) { std::swap(etai, etat); }
	// Compute sini using Snell's law
	float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
	// Total internal reflection
	if (sint >= 1) {
		return 1.f;
	}
	else {
		float cost = sqrtf(std::max(0.f, 1 - sint * sint));
		cosi = fabsf(cosi);
		float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
		float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
		return (Rs * Rs + Rp * Rp) / 2;
	}
	// As a consequence of the conservation of energy, transmittance is given by:
	// kt = 1 - kr;
}
glm::vec3 refract(const glm::vec3& incident, const glm::vec3& normal, const float& ior)
{
	float cosi = glm::clamp(glm::dot(incident, normal), -1.f, 1.f);
	float etai = 1, etat = ior;
	glm::vec3 n = normal;
	if (cosi < 0) { cosi = -cosi; }
	else { std::swap(etai, etat); n = -normal; }
	float eta = etai / etat;
	float k = 1 - eta * eta * (1 - cosi * cosi);
	return k < 0 ? glm::vec3(0.f) : eta * incident + (eta * cosi - sqrtf(k)) * n;
}
glm::vec3 reflect(const glm::vec3& incident, const glm::vec3& normal)
{
	return incident - 2 * glm::dot(incident, normal) * normal;
}
glm::vec3 Whitted::TraceRay(Ray &ray, int depth)
{
	if (depth > maxDepth)
	{
		return glm::vec3(0.4f);
	}
	glm::vec3 pointLight(1.f, 1.f, 1.f);
	float intensity = 10.f;
	
	float closest = std::numeric_limits<float>::infinity();
	SurfaceInteraction intersection;
	bool hit = m_ActiveScene->Intersect(ray, closest, intersection);

	if (!hit)
	{
		return glm::vec3(0.4f);
	}
	else
	{
		glm::vec3 surfaceNormal = intersection.hit_normal;
		glm::vec3 hitPoint = ray(closest);
		glm::vec3 hitColor = m_ActiveScene->materials[intersection.material].albedo;
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
				float pathLength = glm::dot(path, path);
				glm::vec3 attentuation = glm::exp(-(glm::vec3(1.f)-hitColor)* .99f *pathLength);
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
		else if (m_ActiveScene->materials[intersection.material].metallic)
		{
			glm::vec3 refl = glm::normalize(reflect(ray.direction, surfaceNormal));
			glm::vec3 reflO = glm::dot(refl, surfaceNormal) < 0 ? hitPoint - surfaceNormal * 0.0001f : hitPoint + surfaceNormal * 0.0001f;
			Ray reflected(reflO, refl);
			return hitColor*TraceRay(reflected, depth + 1);
		}
		else
		{
			glm::vec3 directionToLight = pointLight - hitPoint;
			float dist2 = glm::dot(directionToLight, directionToLight);
			
			directionToLight = glm::normalize(directionToLight);
			Ray shadow(hitPoint+ surfaceNormal *0.0001f, directionToLight);
			float occluded = std::numeric_limits<float>::infinity();
			SurfaceInteraction nop;
			if (!m_ActiveScene->Intersect(shadow, occluded, nop))
			{
				float dist = glm::sqrt(dist2);
				if (occluded < dist) {
					return glm::vec3(0.f);
				}
			}
			float diffuse = intensity * dot(surfaceNormal, directionToLight) * glm::one_over_pi<float>() / dist2;
			return glm::clamp(diffuse * hitColor, 0.f, 1.f);
			
			/*else
				return  glm::vec3(0.f);*/
			// ambient lighting based on sky color
			/*Ray skylight(hitPoint + surfaceNormal * 0.00001f, surfaceNormal);
			if (m_ActiveScene->Intersect(skylight, occluded) < 0)
			{
				color += hitColor * 0.4f * glm::one_over_pi<float>();
			}
			return color;*/
		}
		/*glm::vec3 albedo = m_ActiveScene->spheres[objIndex].GetAlbedo(ray(closest));
		glm::vec3 normal = m_ActiveScene->spheres[objIndex].GetNormal(ray(closest));
		glm::vec3 directionToLight = glm::normalize(pointLight - ray(closest));
		float cosine = glm::dot(directionToLight, normal)*glm::one_over_pi<float>();
		glm::vec4 color(intensity * albedo * cosine, 1.f);
		color = glm::clamp(color, 0.f, 1.f);

		return color;*/
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