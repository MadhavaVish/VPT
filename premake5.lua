-- premake5.lua
workspace "VPT"
   architecture "x64"
   configurations { "Debug", "Release", "Dist" }
   startproject "VPT"

outputdir = "%{cfg.buildcfg}-%{cfg.system}-%{cfg.architecture}"
include "Walnut/WalnutExternal.lua"

include "VPT"