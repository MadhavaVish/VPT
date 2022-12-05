#include "Walnut/Application.h"
#include "Walnut/EntryPoint.h"

#include "Walnut/Image.h"
#include "Walnut/Timer.h"
#define GLM_FORCE_DEFAULT_ALIGNED_GENTYPES
#define GLM_FORCE_AVX

#include "Integrators/WhittedRenderer.hpp"
#include "Integrators/PathTracer.hpp"
#include "Utils/Camera.hpp"
#include <memory>
using namespace Walnut;

class ExampleLayer : public Walnut::Layer
{
public:
	ExampleLayer() : m_Camera(45.0f, 0.1f, 100.0f) {}

	virtual void OnUpdate(float ts) override {
		if (m_Camera.OnUpdate(ts))
		{
			m_Renderer.Reset();
		}
	}
	virtual void OnUIRender() override
	{
		ImGui::Begin("Settings");
		ImGui::Text("Last render time: %.3fms", m_lastRenderTime);
		if (ImGui::Button("Render")) 
		{
			Render();
		}
		ImGui::Checkbox("Accumulate", &m_Renderer.GetSettings().Accumulate);
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
		m_Camera.OnResize(m_ViewportWidth, m_ViewportHeight);
		m_Renderer.Render(m_Scene, m_Camera);
		
		m_lastRenderTime = timer.ElapsedMillis();

	}
private:
	Scene m_Scene;
	//Whitted m_Renderer;
	PathTracer m_Renderer;
	Camera m_Camera;
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