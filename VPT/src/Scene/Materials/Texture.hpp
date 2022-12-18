#pragma once
#include <stb_image.h>
#include <string>
#include <iostream>
#include <glm/glm.hpp>
#include <glm/gtx/common.hpp>
class Texture
{
public:
	const float rgbScale = 1 / 255.f;
	const static int rgb = 3;
	Texture() : data(nullptr), width(0), height(0), bytesPerLine(0) { };
	Texture(const std::string& filepath)
	{
		int channels_per_pixel = rgb;
		data = stbi_load(filepath.c_str(),&width, &height, &channels_per_pixel, channels_per_pixel);

		if (!data)
		{
			std::cerr << "pain lmao" << std::endl;
		}
		bytesPerLine = rgb * width;
	}
	~Texture()
	{
		stbi_image_free(data);
	}
	glm::vec3 sampleImageTexture(const float& u, const float& v) const
	{
		if (u != u || v != v)
		{
			return glm::vec3(0.f);
		}
		if (!data) return glm::vec3(0.5f, 0.1f, 1.f);

		float a = glm::clamp(u, 0.f, 1.f);
		float b =  1 - glm::clamp(v, 0.f, 1.f);

		int i = static_cast<int>(a * width);
		int j = static_cast<int>(b * height);

		if (i >= width) i = width - 1;
		if (j >= height) j = height - 1;

		auto pixel = data + j * bytesPerLine + i * rgb;


		return glm::vec3(rgbScale * pixel[0], rgbScale * pixel[1], rgbScale * pixel[2]);

	}
private:
	unsigned char* data;
	int width, height;
	int bytesPerLine;
};