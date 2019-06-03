cmake_minimum_required(VERSION 3.8)
# The C is verquired to build with VTK when it has MPI
project(ttk VERSION 0.9.8 LANGUAGES CXX C)

option(TTK_BUILD_VTK_WRAPPERS "Build the TTK VTK Wrappers" ON)
option(TTK_BUILD_PARAVIEW_PLUGINS "Build the TTK ParaView Plugins" ON)
option(TTK_BUILD_STANDALONE_APPS "Build the TTK Standalone Applications" ON)

if(NOT TTK_BUILD_VTK_WRAPPERS)
  set(TTK_BUILD_STANDALONE_APPS OFF)
endif()

if(TTK_BUILD_PARAVIEW_PLUGINS)
  set(TTK_BUILD_VTK_WRAPPERS ON)
endif()

# find dependencies
include(config.cmake)

# Functiond to read ttk.module files
if(VTK_FOUND)
  # required for both wrappers and paraview
  include(CMake/VTKModule.cmake)
endif()

# Base code
# ---------

include(CMake/BaseCode.cmake)
add_subdirectory(core)

# Paraview plugin
# ---------------

if(TTK_BUILD_PARAVIEW_PLUGINS)
  # This option allows to build the TopologyToolKit.so file for paraview
  option(BUILD_SHARED_LIBS "Build TTK as shared lib" ON)
  include(CMake/ParaViewFilter.cmake)
  add_subdirectory(paraview)
endif()

# VTK Wrappers
# ------------

if(TTK_BUILD_VTK_WRAPPERS)
  include(CMake/VTKWrapper.cmake)
  add_subdirectory(core/vtk)
endif()

# Standalones
# -----------

if(TTK_BUILD_STANDALONE_APPS)
  add_subdirectory(standalone)
endif()

# Status
# ------

include(CMake/Print.cmake)
ttk_print_summary()
