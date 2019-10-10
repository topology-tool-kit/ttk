cmake_minimum_required(VERSION 3.12)
# The C is verquired to build with VTK when it has MPI
project(ttk VERSION 0.9.8 LANGUAGES CXX C)

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(FATAL_ERROR "Build in sources is not supported by TTK, please use a separate build folder")
endif()


option(TTK_BUILD_VTK_WRAPPERS "Build the TTK VTK Wrappers" ON)
option(TTK_BUILD_PARAVIEW_PLUGINS "Build the TTK ParaView Plugins" ON)
option(TTK_BUILD_STANDALONE_APPS "Build the TTK Standalone Applications" ON)
option(TTK_WHITELIST_MODE "Explicitely enable each filter" OFF)

# This option allows library to be built dynamic
# like the TopologyToolKit.so file for paraview
option(BUILD_SHARED_LIBS "Build TTK as shared lib" ON)

mark_as_advanced(TTK_WHITELIST_MODE BUILD_SHARED_LIBS)

if(TTK_BUILD_STANDALONE_APPS AND NOT TTK_BUILD_VTK_WRAPPERS)
  message(WARNING "Can't build standalones without the VTK wrappers: disable")
  set(TTK_BUILD_STANDALONE_APPS OFF CACHE BOOL "Build the cmd and gui commands" FORCE)
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
  if(${TTK_WHITELIST_MODE})
    set(TTK_ENABLE_FILTER_DEFAULT OFF CACHE INTERNAL "Default value for each filter build")
  else()
    set(TTK_ENABLE_FILTER_DEFAULT ON CACHE INTERNAL "Default value for each filter build")
  endif()
  include(CMake/ParaViewFilter.cmake)
  add_subdirectory(paraview)
endif()

# VTK Wrappers
# ------------

if(TTK_BUILD_VTK_WRAPPERS)
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
