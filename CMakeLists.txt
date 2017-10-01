cmake_minimum_required(VERSION 3.2)
project(ttk VERSION 0.9.2 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)

option(TTK_ENABLE_KAMIKAZE "Enable Kamikaze compilation mode" FALSE)
option(TTK_ENABLE_CPU_OPTIMIZATION "Enable native CPU optimizations" TRUE)
option(TTK_BUILD_VTK_WRAPPERS "Build the TTK VTK Wrappers" ON)
option(TTK_BUILD_PARAVIEW_PLUGINS "Build the TTK ParaView Plugins" ON)
option(TTK_BUILD_STANDALONE_APPS "Build the TTK Standalone Applications" ON)

add_subdirectory(core)
if (TTK_BUILD_PARAVIEW_PLUGINS)
	add_subdirectory(paraview)
endif()

if (TTK_BUILD_STANDALONE_APPS)
	add_subdirectory(standalone)
endif()

