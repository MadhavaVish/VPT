#include "Walnut/Application.h"
#include "Walnut/EntryPoint.h"

#include "Walnut/Image.h"
#include "Walnut/Timer.h"
#define GLM_FORCE_DEFAULT_ALIGNED_GENTYPES
#define GLM_FORCE_AVX

//#include "Integrators/WhittedRenderer.hpp"
#include "Integrators/PathTracer.hpp"
#include "Integrators/WhittedRenderer.hpp"
#include "Utils/Camera.hpp"
#include "Utils/Camera.hpp"
#include <memory>
using namespace Walnut;

class ExampleLayer : public Walnut::Layer
{
public:
	ExampleLayer() {}

	virtual void OnUpdate(float ts) override {
		if (cam.OnUpdate(ts))
		{
			m_Renderer.Reset();
		}
	}
	virtual void OnUIRender() override
	{
		ImGui::Begin("Settings");
		if (ImGui::CollapsingHeader("Camera Settings"))
		{
			if (ImGui::SliderFloat("Focal Length", &cam.settings.focal_length, 10.f, 200.f))
			{
				cam.OnResize(m_ViewportWidth, m_ViewportHeight);
				m_Renderer.Reset();
			}
			if (ImGui::SliderFloat("F-Stop", &cam.settings.f_stop, 0.01f, 32.f))
			{
				cam.OnResize(m_ViewportWidth, m_ViewportHeight);
				m_Renderer.Reset();
			}
			if (ImGui::SliderFloat("Focal Distance", &cam.settings.focus_dist, 0.1f, 200.f))
			{
				cam.OnResize(m_ViewportWidth, m_ViewportHeight);
				m_Renderer.Reset();
			}
			if (ImGui::SliderFloat("Sensor Width", &cam.settings.sensor_width, 1.f, 36.f))
			{
				cam.OnResize(m_ViewportWidth, m_ViewportHeight);
				m_Renderer.Reset();
			}
		}
		//ImGui::EndTabBar();
		ImGui::Text("Last render time: %.3fms", m_lastRenderTime);
		if (ImGui::Button("Render")) 
		{
			Render();
		}
		
		ImGui::Checkbox("Accumulate", &m_Renderer.GetSettings().Accumulate);
		//ImGui::Checkbox("Vignette", &m_Renderer.GetSettings().Vignette);
		//ImGui::SliderFloat("Amount", &m_Renderer.GetSettings().VignetteAmount, 0.f, 20.f);
		if (ImGui::Button("Reset"))
		{
			m_Renderer.Reset();
		}
		ImGui::End();

		ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0f, 0.0f));
		ImGui::Begin("Viewport");
		m_ViewportWidth = ImGui::GetContentRegionAvail().x;
		m_ViewportHeight = ImGui::GetContentRegionAvail().y;
		auto image = m_Renderer.GetFinalImage();
		if(image)
			ImGui::Image(image->GetDescriptorSet(), { (float)image->GetWidth(), (float)image->GetHeight() }, ImVec2(0,1), ImVec2(1, 0));
		ImGui::PopStyleVar();
		ImGui::End();

		Render();
	}
	void Render()
	{
		Timer timer;
		m_Renderer.OnResize(m_ViewportWidth, m_ViewportHeight);
		cam.OnResize(m_ViewportWidth, m_ViewportHeight);
		m_Renderer.Render(m_Scene, cam);
		
		m_lastRenderTime = timer.ElapsedMillis();

	}
private:
	Scene m_Scene;
	Camera cam;
	//Whitted m_Renderer;
	PathTracer m_Renderer;
	uint32_t  m_ViewportHeight = 0, m_ViewportWidth = 0;
	float m_lastRenderTime = 0.f;
};

Walnut::Application* Walnut::CreateApplication(int argc, char** argv)
{
	Walnut::ApplicationSpecification spec;
	spec.Name = "Walnut Example";

	Walnut::Application* app = new Walnut::Application(spec);
	app->PushLayer<ExampleLayer>();
	app->SetMenubarCallback([app]()
	{
		if (ImGui::BeginMenu("File"))
		{
			if (ImGui::MenuItem("Exit"))
			{
				app->Close();
			}
			ImGui::EndMenu();
		}
	});
	return app;
}