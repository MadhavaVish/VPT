# CPU Path Tracer

This project was built for my Advanced Graphics class. 

## Features
- Unbiased path tracing
- Model loading
- high performance bvh
- TODO top level acceleration structure
- Texturing
- Disney BRDF
- Realistic Camera Model with controls that map to real world cameras

## Getting Started

- Clone recursively to make sure all the dependencies are met. Install vulkan using the sdk or other method as long as premake can find it you're good!
- Go to scripts and run Setup.bat
- Launch the visual studio project and compile!

## Gallery
![image](https://github.com/MadhavaVish/VPT/assets/19480221/9e72477d-7e12-44b7-a723-fdcc31511319)
*Disney BRDF and rough metallics and dielectrics in action! (Current state of the renderer)*
![image](https://github.com/MadhavaVish/VPT/assets/19480221/6573f0bf-505b-41c9-8bf5-797d3dc63de7)
*Old and incorrect implementation of defocus blur but still looks cool i think*
![image](https://github.com/MadhavaVish/VPT/assets/19480221/0f3466a7-3d53-4f78-98a2-ba840f21bbbb)
*Emissive globe inside transparent glass ball (not quite a layered material model yet but still fun to mess around with*
![image](https://github.com/MadhavaVish/VPT/assets/19480221/b47a625d-05f5-42e1-80c9-d4cea641f7cf)
*Cornell meet Vincent*
![image](https://github.com/MadhavaVish/VPT/assets/19480221/3907bfb7-7927-44d1-a9ed-541da0773dc6)
*Bread on a desk, this time with a correct DOF implementation and textures!*
![image](https://github.com/MadhavaVish/VPT/assets/19480221/d1cdb834-4c53-42f8-9a5b-e26001951820)
*Using Jacco Bikkers articles a BVH handling a 28 million triangle Lucy!*
![image](https://github.com/MadhavaVish/VPT/assets/19480221/ce0dd3ec-6093-409b-b477-8b1b42338714)
*Even without translucency or multiple material maps it turned out quite pretty*

