-- premake5.lua
workspace "VPT"
   architecture "x64"
   configurations { "Debug", "Release", "Dist" }
   startproject "VPT"

includedirs
{
   "Walnut/vendor/stb_image"
}
openmp "On"
outputdir = "%{cfg.buildcfg}-%{cfg.system}-%{cfg.architecture}"
include "Walnut/WalnutExternal.lua"

include "VPT"