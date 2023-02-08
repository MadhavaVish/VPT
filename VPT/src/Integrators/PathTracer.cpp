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

	//m_ActiveCamera = &camera;
	//std::for_each(std::execution::par, m_VerticalIter.begin(), m_VerticalIter.end(),
	//	[this](uint32_t y) {
	//		for (int x = 0; x < m_FinalImage->GetWidth(); x++)
	//		{
	//			//m_ActiveCamera = &camera;
	//			Ray ray;
	//			glm::vec2 randomOffset(Walnut::Random::Float(), Walnut::Random::Float());
	//			ray = m_ActiveCamera->getPrimaryRay(x, y, randomOffset);

	//			glm::vec3 color(TraceRay(ray,0));
	//			accumulator[x + y * m_FinalImage->GetWidth()] += color;

	//			glm::vec4 accumulated(glm::sqrt(accumulator[x + y * m_FinalImage->GetWidth()] / (float)frameIndex), 1.f);
	//			accumulated = glm::clamp(accumulated, glm::vec4(0.f), glm::vec4(1.f));

	//			m_ImageData[x + y * m_FinalImage->GetWidth()] = Utils::ConvertToRGBA(accumulated);
	//		}
	//	});
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

			glm::vec4 accumulated(glm::pow(accumulator[x + y * m_FinalImage->GetWidth()] / (float)frameIndex, glm::vec3(1.f/2.2f)), 1.f);
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
static glm::vec3 lerp(glm::vec3 x, glm::vec3 y, float t) {
	return x * (1.f - t) + y * t;
}
static float max3(glm::vec3 v) {
	return glm::max(glm::max(v.x, v.y), v.z);
}
static float computeRelativeIOR(const Material& mat, const glm::vec3& wo) {
	bool entering = cosTheta(wo) > 0.0f;
	return entering ? 1.0f / mat.ior_ : mat.ior_;
}
static void computeLobeProbabilities(const Material& mat, const glm::vec3& wo, float& pDiffuse, float& pSpecular, float& pTransmission) {
	float eta = computeRelativeIOR(mat, wo);
	glm::vec3 f0 = lerp(schlickF0FromRelativeIOR(eta), mat.baseColor_, mat.metalness_);
	glm::vec3 fresnel = Fr_Schlick(abs(cosTheta(wo)), f0);

	float diffuseWeight = (1.0f - mat.metalness_) * (1.0f - mat.transmission_);
	float transmissionWeight = (1.0f - mat.metalness_) * mat.transmission_;

	pDiffuse = max3(mat.baseColor_) * diffuseWeight;
	pSpecular = max3(fresnel);
	pTransmission = max3(glm::vec3(1.0f) - fresnel) * transmissionWeight;

	float normFactor = 1.0f / (pDiffuse + pSpecular + pTransmission);
	pDiffuse *= normFactor;
	pSpecular *= normFactor;
	pTransmission *= normFactor;
}
static float pdf(const Material& mat, const glm::vec3& wi, const glm::vec3& wo) {
	float eta = computeRelativeIOR(mat, wo);

	float pDiffuse, pSpecular, pTransmission;
	computeLobeProbabilities(mat, wo, pDiffuse, pSpecular, pTransmission);

	return pDiffuse * pdfCosineHemisphere(wi, wo)
		+ pSpecular * pdfGGX_VNDF_reflection(wi, wo, mat.alpha_)
		+ pTransmission * pdfGGX_VNDF_transmission(wi, wo, eta, mat.alpha_);
}


static glm::vec3 sampleDirection(const Material& mat, const glm::vec3& wo, float r, float u1, float u2, float* pdf) {
	float eta = computeRelativeIOR(mat, wo);
	float pDiffuse, pSpecular, pTransmission;
	computeLobeProbabilities(mat, wo, pDiffuse, pSpecular, pTransmission);

	glm::vec3 wi;
	if (r < pDiffuse) {

		wi = sign(cosTheta(wo)) * sampleCosineHemisphere(u1, u2);
		//assert(isNormalized(wi));
	}
	else if (r < pDiffuse + pSpecular) {
		glm::vec3 wo_upper = sign(cosTheta(wo)) * wo; // sign(+wo) * +wo = +wo, sign(-wo) * -wo = +wo
		glm::vec3 wh = sign(cosTheta(wo)) * sampleGGX_VNDF(wo_upper, mat.alpha_, u1, u2);
		if (dot(wo, wh) < 0.0f) {
			return glm::vec3(0.0f);
		}

		wi = reflect(wo, wh);
		if (!sameHemisphere(wi, wo)) {
			wi = -wi;
			//return glm::vec3(0.0f);
		}
	}
	else {

		glm::vec3 wo_upper = sign(cosTheta(wo)) * wo; // sign(+wo) * +wo = +wo, sign(-wo) * -wo = +wo
		glm::vec3 wh = sign(cosTheta(wo)) * sampleGGX_VNDF(wo_upper, mat.alpha_, u1, u2);

		refract(wo, wh, eta, wi, mat.ior_);
	}

	if (pdf) {
		*pdf = pDiffuse * pdfCosineHemisphere(wi, wo)
			+ pSpecular * pdfGGX_VNDF_reflection(wi, wo, mat.alpha_)
			+ pTransmission * pdfGGX_VNDF_transmission(wi, wo, eta, mat.alpha_);
	}

	return wi;
}

static glm::vec3 evaluate(const Material& mat, const glm::vec3 hitColor, const glm::vec3& wi, const glm::vec3& wo) {
	float eta = computeRelativeIOR(mat, wo);
	glm::vec3 f0 = lerp(schlickF0FromRelativeIOR(eta), hitColor, mat.metalness_);

	glm::vec3 diffuse = diffuse_Lambert(wi, wo, hitColor);
	glm::vec3 specular = microfacetReflection_GGX(wi, wo, f0, eta, mat.alpha_);
	glm::vec3 transmission = hitColor * microfacetTransmission_GGX(wi, wo, f0, eta, mat.alpha_);

	float diffuseWeight = (1.0f - mat.metalness_) * (1.0f - mat.transmission_);
	float transmissionWeight = (1.0f - mat.metalness_) * mat.transmission_;

	return diffuseWeight * diffuse + specular + transmissionWeight * transmission;
}

glm::vec3 PathTracer::TraceRay(Ray& ray, int depth)
{
	int bounce = 0;
	glm::vec3 throughput(1.f);
	while(true)
	{
		Intersection isect;
		if (m_ActiveScene->Intersect(ray, isect))
		{
			SurfaceInteraction interaction = m_ActiveScene->getSurfaceProperties(ray, isect);
			Frame surf(interaction.hit_normal);

			Material mat = m_ActiveScene->materials[interaction.materialIdx];
			glm::vec3 hitColor = mat.baseColor_;
			float rrProb = 1;
			if (bounce > 3)
			{
				rrProb = glm::clamp(std::max({ throughput.x, throughput.y, throughput.z }), 0.05f, 0.95f);
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
			if (glm::dot(mat.emittance_,mat.emittance_) > 0.f)
			{
				return throughput * hitColor * mat.emittance_;
			}
			float pdf;
			glm::vec3 dir = sampleDirection(mat, surf.ToLocal(-ray.direction) , Walnut::Random::Float(), Walnut::Random::Float(), Walnut::Random::Float(), &pdf);
			float cosThetaI = abs(cosTheta(dir));
			if (cosThetaI <= 0.0f || pdf <= 0.0f) {
				break;
			}
			glm::vec3 bsdf = evaluate(mat, hitColor, dir, surf.ToLocal(-ray.direction));
			ray = getRay(rayPnt(ray, isect.t_hit) + sign(cosTheta(dir)) * interaction.hit_normal * 0.001f, surf.ToWorld(dir));
			throughput *= ((bsdf * cosThetaI)/(pdf*rrProb));
			bounce++;
			continue;
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

