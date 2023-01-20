#pragma once
#include <glm/glm.hpp>
#include "../../Utils/Ray.hpp"
#include "../../Utils/Utils.hpp"
#include <Walnut/Random.h>
#include <glm/gtc/constants.hpp>

struct Material {
	glm::vec3 albedo;
	float radiance;
	float exponent;
	bool mirror;
	bool glass;
	float ks;
	float ior;
	int textureIndex = -1;
};

static glm::vec3 sampleCosineHemisphere(float& pdf)
{
	float a = glm::two_pi<float>() * Walnut::Random::Float();
	float b = Walnut::Random::Float();
	float c = glm::sqrt(1 - b);
	float d = glm::sqrt(b);
	pdf = d * glm::one_over_pi<float>();
	return glm::vec3(glm::cos(a) * c, glm::sin(a) * c, d);
}
static float fresnel(const glm::vec3& I, const glm::vec3& N, const float& ior)
{
	float cosi = glm::clamp(glm::dot(I, N), -1.f, 1.f);
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
static glm::vec3 refract(const glm::vec3& incident, const glm::vec3& normal, const float& ior)
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
static glm::vec3 reflect(const glm::vec3& incident, const glm::vec3& normal)
{
	return incident - 2 * glm::dot(incident, normal) * normal;
}

static glm::vec3 SampleUniformHemisphere(float& pdf)
{
	float a = Walnut::Random::Float();
	float b = Walnut::Random::Float();
	float cosTheta = a;
	float sinTheta = glm::sqrt(1 - cosTheta * cosTheta);
	float cosPhi = glm::cos(glm::two_pi<float>() * b);
	float sinPhi = glm::sin(glm::two_pi<float>() * b);
	pdf = glm::one_over_two_pi<float>();
	return glm::vec3(cosPhi * sinTheta, sinPhi * sinTheta, cosTheta);
}

static glm::vec3 samplePowerCosineHemisphere(const float& exponent, float& pdf)
{
	float a = glm::two_pi<float>() * Walnut::Random::Float();
	//float theta = glm::acos(glm::pow(1 - Walnut::Random::Float(), 1 / (1 + exponent)));
	float b = glm::pow(Walnut::Random::Float(), 1 / (exponent + 1));
	float c = glm::sqrt(1 - b * c);
	//float x = glm::sin(theta) * glm::cos(a);
	//float y = glm::sin(theta) * glm::sin(a);
	//float z = glm::sqrt(1 - y * x);

	pdf = (exponent + 1) * glm::pow(b, exponent) * glm::one_over_two_pi<float>();
	//pdf = 1;

	return glm::vec3(glm::cos(a) * c, glm::sin(a) * c, b);
	//return glm::vec3(x, y, z);
}

static glm::vec2 SampleConcentricDisc()
{
	float phi, r;

	float a = 2 * Walnut::Random::Float() - 1;   /* (a,b) is now on [-1,1]^2 */
	float b = 2 * Walnut::Random::Float() - 1;

	if (a > -b)      /* region 1 or 2 */
	{
		if (a > b)   /* region 1, also |a| > |b| */
		{
			r = a;
			phi = (glm::pi<float>() / 4.f) * (b / a);
		}
		else        /* region 2, also |b| > |a| */
		{
			r = b;
			phi = (glm::pi<float>() / 4.f) * (2.f - (a / b));
		}
	}
	else            /* region 3 or 4 */
	{
		if (a < b)   /* region 3, also |a| >= |b|, a != 0 */
		{
			r = -a;
			phi = (glm::pi<float>() / 4.f) * (4.f + (b / a));
		}
		else        /* region 4, |b| >= |a|, but a==0 and b==0 could occur. */
		{
			r = -b;

			if (b != 0)
				phi = (glm::pi<float>() / 4.f) * (6.f - (a / b));
			else
				phi = 0;
		}
	}

	glm::vec2 res;
	res.x = r * std::cos(phi);
	res.y = r * std::sin(phi);
	return res;
}