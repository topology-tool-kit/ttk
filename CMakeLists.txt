cmake_minimum_required(VERSION 3.2)
project(ttk VERSION 0.9.4 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)

# Set a predefined build type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
   message(STATUS "Setting build type to 'Release'.")
   set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
   set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

option(TTK_USE_64BIT_IDS "Enable processing on large datasets" OFF)
option(TTK_ENABLE_KAMIKAZE "Enable Kamikaze compilation mode" OFF)
option(TTK_ENABLE_CPU_OPTIMIZATION "Enable native CPU optimizations" ON)

option(TTK_BUILD_VTK_WRAPPERS "Build the TTK VTK Wrappers" ON)
option(TTK_BUILD_PARAVIEW_PLUGINS "Build the TTK ParaView Plugins" ON)
option(TTK_BUILD_STANDALONE_APPS "Build the TTK Standalone Applications" ON)

if (TTK_BUILD_STANDALONE_APPS)
	set(TTK_BUILD_VTK_WRAPPERS ON)
endif()

add_subdirectory(core)
if (TTK_BUILD_PARAVIEW_PLUGINS)
	add_subdirectory(paraview)
endif()

if (TTK_BUILD_STANDALONE_APPS)
	add_subdirectory(standalone)
endif()

ttk_print_summary()
