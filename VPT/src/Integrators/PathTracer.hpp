#pragma once
#include "../Scene/Scene.hpp"

#include "../Utils/Camera.hpp"
#include "../Utils/Ray.hpp"

#include <Walnut/Image.h>
#include <glm/glm.hpp>
#include <memory>

class PathTracer
{
public:
	struct Settings
	{
		bool Accumulate = true;
	};
	PathTracer() = default;

	void Render(const Scene& scene, const Camera& camera);
	void OnResize(uint32_t width, uint32_t height);
	void Reset() { frameIndex = 1; }
	std::shared_ptr<Walnut::Image> GetFinalImage() const { return m_FinalImage; }
	Settings& GetSettings() { return settings; }
private:
	glm::vec3 TraceRay(Ray& ray);
private:
	const Scene* m_ActiveScene;
	const Camera* m_ActiveCamera;
	std::shared_ptr<Walnut::Image> m_FinalImage;
	Settings settings;
	uint32_t* m_ImageData = nullptr;
	glm::vec3* accumulator = nullptr;
	uint32_t frameIndex = 1;
};