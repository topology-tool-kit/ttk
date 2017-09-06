cmake_minimum_required(VERSION 3.5)
project(ttk)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(TTK_BUILD_PARAVIEW_PLUGINS "Build the TTK ParaView Plugins" OFF)
if (NOT TTK_BUILD_PARAVIEW_PLUGINS)
	option(TTK_BUILD_VTK_WRAPPERS "Build the TTK VTK Wrappers" OFF)
	find_package(VTK REQUIRED)
	include(${VTK_USE_FILE})
else()
	set(TTK_BUILD_VTK_WRAPPERS ON CACHE BOOL "Build the TTK VTK Wrappers" FORCE)
	find_package(ParaView REQUIRED)
	include(${PARAVIEW_USE_FILE})
endif()

add_subdirectory(core)
add_subdirectory(standalone)

