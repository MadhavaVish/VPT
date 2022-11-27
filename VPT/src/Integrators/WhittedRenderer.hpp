#pragma once
#include "../Scene/Scene.hpp"

#include "../Utils/Camera.hpp"
#include "../Utils/Ray.hpp"

#include <Walnut/Image.h>
#include <glm/glm.hpp>
#include <memory>
class Whitted {
public:
	Whitted() = default;
	//Whitted(Scene* scene): scene(scene){}
	
	void Render(const Scene &scene, const Camera &camera);
	void OnResize(uint32_t width, uint32_t height);
	std::shared_ptr<Walnut::Image> GetFinalImage() const { return m_FinalImage; }
private:
	glm::vec4 TraceRay(int x, int y);
private:
	const Scene* m_ActiveScene;
	const Camera* m_ActiveCamera;
	std::shared_ptr<Walnut::Image> m_FinalImage;
	uint32_t* m_ImageData = nullptr;
};