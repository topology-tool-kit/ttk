cmake_minimum_required(VERSION 3.2)
project(ttk VERSION 0.9.8 LANGUAGES CXX)

option(TTK_BUILD_VTK_WRAPPERS "Build the TTK VTK Wrappers" ON)
option(TTK_BUILD_PARAVIEW_PLUGINS "Build the TTK ParaView Plugins" ON)
option(TTK_BUILD_STANDALONE_APPS "Build the TTK Standalone Applications" ON)

include(functions.cmake)

if (NOT TTK_BUILD_VTK_WRAPPERS)
  set(TTK_BUILD_STANDALONE_APPS OFF)
endif()

add_subdirectory(core)

if (TTK_BUILD_PARAVIEW_PLUGINS)
  add_subdirectory(paraview)
endif()

if (TTK_BUILD_STANDALONE_APPS)
  add_subdirectory(standalone)
endif()

ttk_print_summary()
