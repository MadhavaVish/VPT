-- premake5.lua
workspace "VPT"
   architecture "x64"
   configurations { "Debug", "Release", "Dist" }
   startproject "VPT"

includedirs
{
   "Walnut/vendor/stb_image",
   "vendor/tiny_obj_loader",
}
openmp "On"
outputdir = "%{cfg.buildcfg}-%{cfg.system}-%{cfg.architecture}"
include "Walnut/WalnutExternal.lua"

include "VPT"