#pragma once
#include "../Scene/Scene.hpp"

#include "../Utils/Camera.hpp"
#include "../Utils/Ray.hpp"

#include <Walnut/Image.h>
#include <glm/glm.hpp>
#include <memory>
class Whitted {
public:
	struct Settings
	{
		bool Accumulate = true;
	};
	Whitted() = default;

	void Render(const Scene &scene, const Camera &camera);
	void OnResize(uint32_t width, uint32_t height);
	void Reset() {};
	Settings& GetSettings() { return settings; }
	std::shared_ptr<Walnut::Image> GetFinalImage() const { return m_FinalImage; }
private:
	glm::vec3 TraceRay(Ray & ray, int depth);
private:
	const Scene* m_ActiveScene;
	const Camera* m_ActiveCamera;
	Settings settings;
	std::shared_ptr<Walnut::Image> m_FinalImage;
	uint32_t* m_ImageData = nullptr;
};