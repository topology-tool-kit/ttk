cmake_minimum_required(VERSION 3.12)
# The C is verquired to build with VTK when it has MPI
project(ttk VERSION 0.9.9 LANGUAGES CXX C)

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(FATAL_ERROR "Build in sources is not supported by TTK, please use a separate build folder")
endif()

option(TTK_BUILD_VTK_WRAPPERS "Build the TTK VTK Wrappers" ON)
option(TTK_BUILD_PARAVIEW_PLUGINS "Build the TTK ParaView Plugins" ON)
option(TTK_BUILD_STANDALONE_APPS "Build the TTK Standalone Applications" ON)
option(TTK_WHITELIST_MODE "Explicitely enable each filter" OFF)
mark_as_advanced(TTK_WHITELIST_MODE BUILD_SHARED_LIBS)

if(${TTK_WHITELIST_MODE})
  set(TTK_ENABLE_FILTER_DEFAULT "DONT_WANT" CACHE INTERNAL "Default value for each filter build")
else()
  set(TTK_ENABLE_FILTER_DEFAULT "WANT" CACHE INTERNAL "Default value for each filter build")
endif()

# This option allows library to be built dynamic
# like the TopologyToolKit.so file for paraview
option(BUILD_SHARED_LIBS "Build TTK as shared lib" ON)

if(TTK_BUILD_STANDALONE_APPS AND NOT TTK_BUILD_VTK_WRAPPERS)
  message(WARNING "Can't build standalones without the VTK wrappers: disable")
  set(TTK_BUILD_STANDALONE_APPS OFF CACHE BOOL "Build the cmd and gui commands" FORCE)
endif()

# find dependencies
include(config.cmake)

# Functiond to read ttk.module files
if(VTK_FOUND)
  # required for both wrappers and paraview
  # Note, may also use it for basecode it we want a base.module file
  include(CMake/VTKModule.cmake)
endif()

# Base code
# ---------

include(CMake/BaseCode.cmake)
add_subdirectory(core)

add_library(baseAll INTERFACE)
add_library(ttkAll INTERFACE)

target_link_libraries(baseAll INTERFACE ${TTK_BASECODES})

install(
  TARGETS
    baseAll
  EXPORT
    TTKBaseTargets
  RUNTIME DESTINATION
    bin/ttk
  ARCHIVE DESTINATION
    lib/ttk
  LIBRARY DESTINATION
    lib/ttk
  )

install(
  EXPORT
    TTKBaseTargets
  NAMESPACE
    ttk::base::
  DESTINATION
    lib/cmake/ttk
  )

configure_file(core/base/TTKBaseConfig.cmake.in TTKBaseConfig.cmake @ONLY)
install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/TTKBaseConfig.cmake"
  DESTINATION
    lib/cmake/ttk
  )

# Paraview plugin
# ---------------

if(TTK_BUILD_PARAVIEW_PLUGINS)
  include(CMake/ParaViewFilter.cmake)

  # Install location
  if(NOT "$ENV{PV_PLUGIN_PATH}" STREQUAL "")
    set(TTK_INSTALL_PLUGIN_DIR
      $ENV{PV_PLUGIN_PATH}
      CACHE
      PATH
        "Directory where the ParaView plugin will be installed"
      )
  else()
    message(WARNING "Please, set the PV_PLUGIN_PATH environement variable, the TTK ParaView plugin will be built but not installed.")
    set(TTK_INSTALL_PLUGIN_DIR
      ""
      CACHE
      PATH
      "Directory where the ParaView plugin will be installed"
      )
  endif()

  add_subdirectory(paraview)
endif()

# VTK Wrappers
# ------------

if(TTK_BUILD_VTK_WRAPPERS)
  add_subdirectory(core/vtk)

  target_link_libraries(ttkAll INTERFACE ${TTK_ENABLED_MODULES})

  # Use folder as configured by the distribution
  include(GNUInstallDirs)
  install(
    TARGETS
      baseAll
      ttkAll
    EXPORT
      TTKVTKTargets
    RUNTIME DESTINATION
      ${CMAKE_INSTALL_BINDIR}/ttk
    ARCHIVE DESTINATION
      ${CMAKE_INSTALL_LIBDIR}/ttk
    LIBRARY DESTINATION
      ${CMAKE_INSTALL_LIBDIR}/ttk
    )

  install(
    EXPORT
      TTKVTKTargets
    NAMESPACE
      ttk::vtk::
    DESTINATION
      lib/cmake/ttk
    )

  configure_file(core/vtk/TTKVTKConfig.cmake.in TTKVTKConfig.cmake @ONLY)
  install(
    FILES
      "${CMAKE_CURRENT_BINARY_DIR}/TTKVTKConfig.cmake"
    DESTINATION
      lib/cmake/ttk
    )
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
