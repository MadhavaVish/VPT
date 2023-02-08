#pragma once
#include <glm/glm.hpp>
#include "../../Utils/Ray.hpp"
#include "../../Utils/Utils.hpp"
#include <Walnut/Random.h>
#include <glm/gtc/constants.hpp>

struct Material {
	/*glm::vec3 albedo;
	float radiance;
	float exponent;
	bool mirror;
	bool glass;
	float ks;
	float ior;*/
	int textureIndex = -1;
    glm::vec3 baseColor_;
    float roughness_;
    float metalness_;
    float transmission_;
    float ior_;
    glm::vec3 emittance_;

    // Precomputed values
    float alpha_;
};

static glm::vec3 sampleCosineHemisphere(float& pdf)
{
	float a = glm::two_pi<float>() * Walnut::Random::Float();
	float b = Walnut::Random::Float();
	float c = glm::sqrt(1 - b);
	float d = glm::sqrt(b);
	pdf = d ;
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
//static glm::vec3 refract(const glm::vec3& incident, const glm::vec3& normal, const float& ior)
//{
//	float cosi = glm::clamp(glm::dot(incident, normal), -1.f, 1.f);
//	float etai = 1, etat = ior;
//	glm::vec3 n = normal;
//	if (cosi < 0) { cosi = -cosi; }
//	else { std::swap(etai, etat); n = -normal; }
//	float eta = etai / etat;
//	float k = 1 - eta * eta * (1 - cosi * cosi);
//	return k < 0 ? glm::vec3(0.f) : eta * incident + (eta * cosi - sqrtf(k)) * n;
//}
//static glm::vec3 reflect(const glm::vec3& incident, const glm::vec3& normal)
//{
//	return incident - 2 * glm::dot(incident, normal) * normal;
//}

static inline bool refract(const glm::vec3& wi, const glm::vec3& n, float eta, glm::vec3& wt, float ior) {
    float cosThetaI = dot(n, wi);
    float sin2ThetaI = glm::max(0.f, 1.f - cosThetaI * cosThetaI);
    float sin2ThetaT = eta * eta * sin2ThetaI;
    if (sin2ThetaT >= 1.f) {
        return false;
    }
    float cosThetaT = std::sqrt(1.f - sin2ThetaT);
    wt = eta * -wi + (eta * cosThetaI - cosThetaT) * n;
    return true;
 //   float cosi = glm::clamp(glm::dot(wi, n), -1.f, 1.f);
	//float etai = 1, etat = ior;
	//glm::vec3 no = n;
 //   if (cosi < 0) { cosi = -cosi; }
	//else { std::swap(etai, etat); no = -n; }
	//float etaa = etai / etat;
	//float k = 1 - eta * eta * (1 - cosi * cosi);
 //   wt = etaa * -wi + (eta * cosi - sqrtf(k)) * no;
 //   return true;
	//return k < 0 ? true : false;
}
static glm::vec3 reflect(const glm::vec3& w, const glm::vec3& n) {
    return static_cast<float>(2) * glm::dot(w, n) * n - w;
}
template <typename T>
static constexpr T remap(T value, T low1, T high1, T low2, T high2) {
    T remapped = low2 + (value - low1) * (high2 - low2) / (high1 - low1);
    return glm::clamp(remapped, low2, high2);
}

template <typename T>
static constexpr T oneMinusEpsilon = static_cast<T>(1.0) - std::numeric_limits<T>::epsilon();

template <typename T, std::enable_if_t<std::is_floating_point<T>::value, bool> = true>
static inline T sign(T a) {
    return std::copysign(static_cast<T>(1), a);
}

//template <typename T>

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
	//return glm::glm::vec3(x, y, z);
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
inline float cosTheta(const glm::vec3& w) { return w.z; }
inline float cosTheta2(const glm::vec3& w) { return w.z * w.z; }
inline float sinTheta2(const glm::vec3& w) { return glm::max(0.0f, 1.0f - cosTheta2(w)); }
inline float sinTheta(const glm::vec3& w) { return std::sqrt(sinTheta2(w)); }
inline float tanTheta(const glm::vec3& w) { return sinTheta(w) / cosTheta(w); }
inline float tanTheta2(const glm::vec3& w) { return sinTheta2(w) / cosTheta2(w); }

inline bool sameHemisphere(const glm::vec3& w, const glm::vec3& wp) {
    return w.z * wp.z > 0.0f;
}


inline glm::vec3 schlickF0FromRelativeIOR(float eta) {
    float a = (1 - eta) / (1 + eta);
    return glm::vec3(a * a);
}

inline glm::vec3 Fr_Schlick(float cosThetaI, const glm::vec3& f0) {
    float a = glm::max(0.0f, 1.0f - cosThetaI);
    float a2 = a * a;
    float a5 = a2 * a2 * a;
    return f0 + (glm::vec3(1.0f) - f0) * a5;
}

inline float D_GGX(const glm::vec3& wh, float alpha) {
    float alpha2 = alpha * alpha;
    float a = 1.0f + cosTheta2(wh) * (alpha2 - 1.0f);
    return alpha2 / (glm::pi<float>() *a * a);
}

inline float G1_Smith_GGX(const glm::vec3& w, float alpha) {
    float tan2ThetaW = tanTheta2(w);
    if (std::isinf(tan2ThetaW)) return 0.0f;
    float alpha2 = alpha * alpha;
    assert(alpha2 * tan2ThetaW >= -1.0f);
    float lambda = (-1.0f + std::sqrt(alpha2 * tan2ThetaW + 1.0f)) / 2.0f;
    return 1.0f / (1.0f + lambda);
}

inline float G2_SmithUncorrelated_GGX(const glm::vec3& wi, const glm::vec3& wo, float alpha) {
    return G1_Smith_GGX(wi, alpha) * G1_Smith_GGX(wo, alpha);
}

inline float G2_SmithHeightCorrelated_GGX(const glm::vec3& wi, const glm::vec3& wo, float alpha) {
    float tan2ThetaO = tanTheta2(wo);
    float tan2ThetaI = tanTheta2(wi);
    if (std::isinf(tan2ThetaO)) return 0.0f;
    if (std::isinf(tan2ThetaI)) return 0.0f;
    float alpha2 = alpha * alpha;
    assert(alpha2 * tan2ThetaO >= -1.0f);
    assert(alpha2 * tan2ThetaI >= -1.0f);
    float lambda_wo = (-1.0f + std::sqrt(alpha2 * tan2ThetaO + 1.0f)) / 2.0f;
    float lambda_wi = (-1.0f + std::sqrt(alpha2 * tan2ThetaI + 1.0f)) / 2.0f;
    return 1.0f / (1.0f + lambda_wo + lambda_wi);
}

inline float G2_None(const glm::vec3& wi, const glm::vec3& wo, float alpha) {
    return 1.0f;
}


// BxDF functions
inline glm::vec3 diffuse_Lambert(const glm::vec3& wi, const glm::vec3& wo, const glm::vec3& diffuseColor) {
    if (!sameHemisphere(wi, wo)) {
        return glm::vec3(0.0f);
    }

    return diffuseColor / glm::pi<float>();
}

inline glm::vec3 microfacetReflection_GGX(const glm::vec3& wi, const glm::vec3& wo, const glm::vec3& f0, float eta, float alpha) {
    if (!sameHemisphere(wi, wo) || cosTheta(wi) == 0.0f || cosTheta(wo) == 0.0f) {
        return glm::vec3(0.0f);
    }

    glm::vec3 wh = wi + wo;
    if (wh.x == 0.0f && wh.y == 0.0f && wh.z == 0.0f) {
        return glm::vec3(0.0f);
    }
    wh = normalize(wh);

    glm::vec3 F;
    if (eta < 1.0f) {
        float cosThetaT = dot(wi, wh);
        float cos2ThetaT = cosThetaT * cosThetaT;
        F = cos2ThetaT > 0.0f ? Fr_Schlick(abs(cosThetaT), f0) : glm::vec3(1.0f);
    }
    else {
        F = Fr_Schlick(abs(dot(wh, wo)), f0);
    }

    float G = G2_SmithHeightCorrelated_GGX(wi, wo, alpha);
    float D = D_GGX(wh, alpha);
    return F * G * D / (4.0f * abs(cosTheta(wi)) * abs(cosTheta(wo)));
}

inline glm::vec3 microfacetTransmission_GGX(const glm::vec3& wi, const glm::vec3& wo, const glm::vec3& f0, float eta, float alpha) {
    if (sameHemisphere(wi, wo) || cosTheta(wi) == 0.0f || cosTheta(wo) == 0.0f) {
        return glm::vec3(0.0f);
    }

    glm::vec3 wh = normalize(wi + eta * wo);
    if (cosTheta(wh) < 0.0f) {
        wh = -wh;
    }

    bool sameSide = dot(wo, wh) * dot(wi, wh) > 0.0f;
    if (sameSide) {
        return glm::vec3(0.0f);
    }

    glm::vec3 F;
    if (eta < 1.0f) {
        float cosThetaT = dot(wi, wh);
        float cos2ThetaT = cosThetaT * cosThetaT;
        F = cos2ThetaT > 0.0f ? Fr_Schlick(abs(cosThetaT), f0) : glm::vec3(1.0f);
    }
    else {
        F = Fr_Schlick(abs(dot(wh, wo)), f0);
    }

    float G = G2_SmithHeightCorrelated_GGX(wi, wo, alpha);
    float D = D_GGX(wh, alpha);
    float denomSqrt = dot(wi, wh) + eta * dot(wo, wh);
    return (glm::vec3(1.0f) - F) * D * G * abs(dot(wi, wh)) * abs(dot(wo, wh))
        / (denomSqrt * denomSqrt * abs(cosTheta(wi)) * abs(cosTheta(wo)));
}


inline glm::vec3 sampleUniformSphere(float u1, float u2) {
    float cosTheta = 1.0f - 2.0f * u1;
    float sinTheta = std::sqrt(glm::max(0.0f, 1.0f - cosTheta * cosTheta));
    float phi = 2.0f * glm::pi<float>() *u2;
    
    return glm::vec3(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
}

inline float pdfUniformSphere(const glm::vec3& wi, const glm::vec3& wo) {
    return 1.0f / (4.0f * glm::pi<float>());
}


inline glm::vec3 sampleCosineHemisphere(float u1, float u2) {
    float cosTheta = std::sqrt(glm::max(0.0f, 1.0f - u1));
    float sinTheta = std::sqrt(u1);
    float phi = 2.0f * glm::pi<float>() *u2;
    return glm::vec3(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
}

inline float pdfCosineHemisphere(const glm::vec3& wi, const glm::vec3& wo) {
    return sameHemisphere(wi, wo) ? cosTheta(wi) / glm::pi<float>() : 0.0f;
}


inline glm::vec3 sampleGGX(float alpha, float u1, float u2) {
    float phi = 2.0f * glm::pi<float>() *u1;
    float cosTheta2 = (1.0f - u2) / ((alpha * alpha - 1.0f) * u2 + 1.0f);
    float sinTheta = std::sqrt(glm::max(0.0f, 1.0f - cosTheta2));
    glm::vec3 wh(sinTheta * std::cos(phi), sinTheta * std::sin(phi), std::sqrt(cosTheta2));
    return wh;
}

inline float pdfGGX_reflection(const glm::vec3& wi, const glm::vec3& wo, float alpha) {
    if (!sameHemisphere(wi, wo)) {
        return 0.0f;
    }

    glm::vec3 wh = normalize(wi + wo);
    float pdf_h = D_GGX(wh, alpha) * abs(cosTheta(wh));
    float dwh_dwi = 1.0f / (4.0f * dot(wi, wh));
    return pdf_h * dwh_dwi;
}

inline float pdfGGX_transmission(const glm::vec3& wi, const glm::vec3& wo, float eta, float alpha) {
    if (sameHemisphere(wi, wo)) {
        return 0.0f;
    }

    glm::vec3 wh = normalize(wi + eta * wo);
    bool sameSide = dot(wo, wh) * dot(wi, wh) > 0.0f;
    if (sameSide) return 0.0f;

    float pdf_h = D_GGX(wh, alpha) * abs(cosTheta(wh));
    float sqrtDenom = dot(wi, wh) + eta * dot(wo, wh);
    float dwh_dwi = abs(dot(wi, wh)) / (sqrtDenom * sqrtDenom);
    return pdf_h * dwh_dwi;
}


// See: http://jcgt.org/published/0007/04/01/paper.pdf
inline glm::vec3 sampleGGX_VNDF(const glm::vec3& wo, float alpha, float u1, float u2) {
    // Transform view direction to hemisphere configuration
    //glm::vec3 woHemi = normalize(glm::vec3(alpha * wo.x, alpha * wo.y, wo.z));

    //// Create orthonormal basis
    //float length2 = woHemi.x * woHemi.x + woHemi.y * woHemi.y;
    //glm::vec3 b1 = length2 > 0.0f
    //    ? glm::vec3(-woHemi.y, woHemi.x, 0.0f) * (1.0f / std::sqrt(length2))
    //    : glm::vec3(1.0f, 0.0f, 0.0f);
    //glm::vec3 b2 = cross(woHemi, b1);

    //// Parameterization of projected area
    //float r = std::sqrt(u1);
    //float phi = 2.0f * glm::pi<float>() *u2;
    //float t1 = r * std::cos(phi);
    //float t2 = r * std::sin(phi);
    //float s = 0.5f * (1.0f + woHemi.z);
    //t2 = (1.0f - s) * std::sqrt(1.0f - t1 * t1) + s * t2;

    //// Reprojection onto hemisphere
    //glm::vec3 whHemi = t1 * b1 + t2 * b2 + std::sqrt(glm::max(0.0f, 1.0f - t1 * t1 - t2 * t2)) * woHemi;

    //// Transforming half vector back to ellipsoid configuration
    //return normalize(glm::vec3(alpha * whHemi.x, alpha * whHemi.y, glm::max(0.0f, whHemi.z)));
    glm::vec3 V = glm::normalize(glm::vec3(alpha * wo.x, alpha * wo.y, wo.z));
    // orthonormal basis
    glm::vec3 T1 = (V.z < 0.9999f) ? glm::normalize(glm::cross(V, glm::vec3(0, 0, 1))) : glm::vec3(1, 0, 0);
    glm::vec3 T2 = cross(T1, V);
    // sample point with polar coordinates (r, phi)
    float a = 1.0 / (1.0 + V.z);
    float r = sqrt(u1);
    float phi = (u2 < a) ? u2 / a * glm::pi<float>() : glm::pi<float>() + (u2 - a) / (1.0f - a) * glm::pi<float>();
    float P1 = r * cos(phi);
    float P2 = r * sin(phi) * ((u2 < a) ? 1.0f : V.z);
    // compute normal
    glm::vec3 N = P1 * T1 + P2 * T2 + glm::sqrt(glm::max(0.0f, 1.0f - P1 * P1 - P2 * P2)) * V;
    // unstretch
    N = glm::normalize(glm::vec3(alpha * N.x, alpha * N.y, glm::max(0.0f, N.z)));
    return N;
}

inline float pdfGGX_VNDF_reflection(const glm::vec3& wi, const glm::vec3& wo, float alpha) {
    if (!sameHemisphere(wi, wo)) {
        return 0.0f;
    }

    glm::vec3 wh = normalize(wi + wo);
    float pdf_h = G1_Smith_GGX(wo, alpha) * D_GGX(wh, alpha) * abs(dot(wh, wo)) / abs(cosTheta(wo));
    float dwh_dwi = 1.0f / (4.0f * dot(wi, wh));
    return pdf_h * dwh_dwi;
}

inline float pdfGGX_VNDF_transmission(const glm::vec3& wi, const glm::vec3& wo, float eta, float alpha) {
    /*if (sameHemisphere(wi, wo)) {
        return pdfGGX_VNDF_reflection(-wi, wo, alpha);
    }*/

    glm::vec3 wh = normalize(wi + eta * wo);
    //bool sameSide = dot(wo, wh) * dot(wi, wh) > 0.0f;
    //if (sameSide) return pdfGGX_VNDF_reflection(-wi, wo, alpha);

    float pdf_h = G1_Smith_GGX(wo, alpha) * D_GGX(wh, alpha) * abs(dot(wh, wo)) / abs(cosTheta(wo));
    float sqrtDenom = dot(wi, wh) + eta * dot(wo, wh);
    float dwh_dwi = abs(dot(wi, wh)) / (sqrtDenom * sqrtDenom);
    return pdf_h * dwh_dwi;
}