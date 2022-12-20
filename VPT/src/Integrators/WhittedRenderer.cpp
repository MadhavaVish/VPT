#include "WhittedRenderer.hpp"
#include "../Utils/Utils.hpp"
#include <glm/gtc/constants.hpp>
#include <iostream>

void Whitted::Render(const Scene &scene, const Camera& camera)
{
	m_ActiveScene = &scene;
	float gamma = 1/2.2f;
	const float invWidth = 1.f / m_FinalImage->GetWidth();
	const float invHeight = 1.f / m_FinalImage->GetWidth();
	#pragma omp parallel for schedule(dynamic)
	for (int y = 0; y < m_FinalImage->GetHeight(); y++)
	{
		m_ActiveCamera = &camera;
		for (uint32_t x = 0; x < m_FinalImage->GetWidth(); x++)
		{
			Ray ray;
			glm::vec2 randomOffset(Walnut::Random::Float() * invWidth, Walnut::Random::Float() * invHeight);
			ray = m_ActiveCamera->getPrimaryRay(x, y, randomOffset);
			glm::vec4 color(glm::pow(TraceRay(ray, 0), glm::vec3(gamma)), 1.f);
			m_ImageData[x +y*m_FinalImage->GetWidth()] = Utils::ConvertToRGBA(color);
		}
	}
	m_FinalImage->SetData(m_ImageData);
}

static int maxBounces = 10;
glm::vec3 Whitted::TraceRay(Ray& ray, int depth)
{
	
	if (depth == maxBounces)
	{
		return glm::vec3(0.f);
	}
	glm::vec3 pointLight(1.f, 1.f, 1.f);
	float intensity = 5.f;
	glm::vec3 color(0.f);
	Intersection isect;
	if (m_ActiveScene->Intersect(ray, isect))
	{
		SurfaceInteraction interaction = m_ActiveScene->getSurfaceProperties(ray, isect);
		Frame surf(interaction.hit_normal);
		
		Material mat = m_ActiveScene->materials[interaction.materialIdx];
		glm::vec3 hitColor = mat.albedo;
		glm::vec3 hitPnt = rayPnt(ray, isect.t_hit);
		if (mat.textureIndex >= 0)
		{
			hitColor *=m_ActiveScene->textures[mat.textureIndex]->sampleImageTexture(interaction.uv.x, interaction.uv.y);
		}
		if (mat.radiance > 0.f)
		{
			return hitColor *mat.radiance;
		}
		if (mat.mirror)
		{
			glm::vec3 dir = surf.ToLocal(-ray.direction);
			dir.x *= -1;
			dir.y *= -1;
			dir = surf.ToWorld(dir);
			ray = getRay(rayPnt(ray, isect.t_hit) + interaction.hit_normal * 0.0001f, dir);

			return hitColor * TraceRay(ray, depth + 1);
		}
		else if (mat.glass)
		{
			glm::vec3 refl = glm::normalize(reflect(ray.direction, interaction.hit_normal));
			glm::vec3 refr = glm::normalize(refract(ray.direction, interaction.hit_normal, mat.ior));
			glm::vec3 reflO = glm::dot(refl, interaction.hit_normal) < 0 ? hitPnt - interaction.hit_normal * 0.0001f : hitPnt + interaction.hit_normal * 0.0001f;
			glm::vec3 refrO = glm::dot(refr, interaction.hit_normal) < 0 ? hitPnt - interaction.hit_normal * 0.0001f : hitPnt + interaction.hit_normal * 0.0001f;
			Ray reflected = getRay(reflO, refl);
			Ray refracted = getRay(refrO, refr);
			float fr = fresnel(ray.direction, interaction.hit_normal, mat.ior);
			if (glm::dot(ray.direction, interaction.hit_normal) > 0)
			{
				glm::vec3 path = hitPnt - ray.origin;
				float pathLength = glm::sqrt(glm::dot(path, path));
				glm::vec3 attentuation = glm::exp(-(hitColor) * 0.6f * pathLength);
				glm::vec3 refrColor = hitColor * attentuation * TraceRay(refracted, depth + 1);
				glm::vec3 reflColor = hitColor * TraceRay(reflected, depth + 1);
				return (fr) * reflColor + (1-fr) * refrColor;
			}
			else
			{
				glm::vec3 refrColor = hitColor*TraceRay(refracted, depth + 1);
				glm::vec3 reflColor = hitColor*TraceRay(reflected, depth + 1);
				return fr * reflColor + (1 - fr) * refrColor;
			}
		}
		else
		{
			
			glm::vec3 directionToLight = pointLight - hitPnt;
			float dist2 = glm::dot(directionToLight, directionToLight);
			directionToLight = glm::normalize(directionToLight);
			Ray shadow = getRay(hitPnt + interaction.hit_normal * 0.0001f, directionToLight);
			//Intersection shad;
			/*if (m_ActiveScene->Intersect(shadow, shad))
			{
				if (shad.t_hit * shad.t_hit < dist2)
				{
					return glm::vec3(0.f);
				}
			}*/
			if (m_ActiveScene->OcclusionBVH(shadow, dist2))
			{
				return glm::vec3(0.f);
			}
			float cosL = glm::dot(interaction.hit_normal, directionToLight);
			glm::vec3 diffColor = hitColor * intensity * cosL * glm::one_over_pi<float>()/ dist2;
			return glm::clamp(diffColor, 0.f, 1.f);
		}
	}
	else
	{
		return m_ActiveScene->getSkyColor(ray);
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