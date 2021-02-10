cmake_minimum_required(VERSION 3.14) # see issue #532
# The C is required to build with VTK when it has MPI
project(ttk VERSION 0.9.10 LANGUAGES CXX C)

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(FATAL_ERROR "Build in sources is not supported by TTK, please use a separate build folder")
endif()

include(CMakeDependentOption)

option(TTK_BUILD_VTK_WRAPPERS "Build the TTK VTK Wrappers" ON)
cmake_dependent_option(TTK_BUILD_PARAVIEW_PLUGINS "Build the TTK ParaView Plugins" ON "TTK_BUILD_VTK_WRAPPERS" OFF)
option(TTK_BUILD_STANDALONE_APPS "Build the TTK Standalone Applications" ON)
option(TTK_WHITELIST_MODE "Explicitely enable each filter" OFF)
mark_as_advanced(TTK_WHITELIST_MODE BUILD_SHARED_LIBS)

# This option allows library to be built dynamic
# like the TopologyToolKit.so file for paraview
option(BUILD_SHARED_LIBS "Build TTK as shared lib" ON)

set(CPACK_PACKAGE_NAME "TTK")
set(CPACK_PACKAGE_FILE_NAME "ttk")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "The Topology ToolKit")
set(CPACK_PACKAGE_VERSION_MAJOR ${TTK_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${TTK_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${TTK_VERSION_PATCH})
set(CPACK_GENERATOR "DEB")
set(CPACK_PACKAGE_CONTACT "Julien Tierny <julien.tierny@sorbonne-universite.fr>")
set(CPACK_PACKAGE_VENDOR "CNRS, Sorbonne University and contributors")
set(CPACK_PACKAGE_HOMEPAGE_URL "https://topology-tool-kit.github.io/")
if(NOT APPLE)
  set(CPACK_RESOURCE_FILE_LICENSE ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE)
  set(CPACK_RESOURCE_FILE_README ${CMAKE_CURRENT_SOURCE_DIR}/README.md)
else()
  # macOS needs license & readme files ending with .txt
  configure_file(LICENSE License.txt COPYONLY)
  configure_file(README.md Readme.txt COPYONLY)
  set(CPACK_RESOURCE_FILE_LICENSE ${PROJECT_BINARY_DIR}/License.txt)
  set(CPACK_RESOURCE_FILE_README ${PROJECT_BINARY_DIR}/Readme.txt)
endif()
set(CPACK_DEBIAN_PACKAGE_DEPENDS
  "ttk-paraview (= 5.8.1), python3-sklearn, libboost-system-dev, python3-dev, libgraphviz-dev, libsqlite3-dev")
# autogenerate dependency information
set (CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
# package will be installed under %ProgramFiles%\${CPACK_PACKAGE_INSTALL_DIRECTORY} on Windows
set(CPACK_PACKAGE_INSTALL_DIRECTORY "TTK")
# let the installer uninstall previous installations on Windows
set(CPACK_NSIS_ENABLE_UNINSTALL_BEFORE_INSTALL ON)
# generate components, fix productbuild packaging for macOS
if(APPLE)
  set(CPACK_COMPONENTS_ALL Unspecified python development)
endif()
# embed Visual C++ redistribuable for Windows + set environment variables
if(WIN32 AND EXISTS "${CMAKE_SOURCE_DIR}/vc_redist.x64.exe")
  install(FILES "${CMAKE_SOURCE_DIR}/scripts/set_windows_env_vars.ps1" DESTINATION bin)
  install(PROGRAMS "${CMAKE_SOURCE_DIR}/vc_redist.x64.exe" DESTINATION bin)
  set(CPACK_VERBATIM_VARIABLES TRUE)
  list(APPEND CPACK_NSIS_EXTRA_INSTALL_COMMANDS
    " ExecWait 'powershell -ExecutionPolicy Bypass -File \"$INSTDIR\\bin\\set_windows_env_vars.ps1\" install'
      ExecWait '$INSTDIR\\bin\\vc_redist.x64.exe /passive /norestart'
      Delete '$INSTDIR\\bin\\vc_redist.x64.exe' ")
  list(APPEND CPACK_NSIS_EXTRA_UNINSTALL_COMMANDS
    " ExecWait 'powershell -ExecutionPolicy Bypass -File \"$INSTDIR\\bin\\set_windows_env_vars.ps1\"' ")
endif()
include(CPack)

if(TTK_BUILD_STANDALONE_APPS AND NOT TTK_BUILD_VTK_WRAPPERS)
  message(WARNING "Can't build standalones without the VTK wrappers: disable")
  set(TTK_BUILD_STANDALONE_APPS OFF CACHE BOOL "Build the cmd and gui commands" FORCE)
endif()

# find dependencies
include(config.cmake)

# Use folder as configured by the distribution
include(GNUInstallDirs)

# Base code
# ---------

include(CMake/CompilerFlags.cmake)
include(CMake/BaseCode.cmake)
add_library(baseAll INTERFACE)

add_subdirectory(core)

target_link_libraries(baseAll INTERFACE ${TTK_BASECODES})

install(
  TARGETS
    baseAll
  EXPORT
    TTKBaseTargets
  RUNTIME DESTINATION
    ${CMAKE_INSTALL_BINDIR}/ttk
  ARCHIVE DESTINATION
    ${CMAKE_INSTALL_LIBDIR}/ttk
  LIBRARY DESTINATION
    ${CMAKE_INSTALL_LIBDIR}/ttk
  )

install(
  EXPORT
    TTKBaseTargets
  NAMESPACE
    ttk::base::
  DESTINATION
    "${CMAKE_INSTALL_LIBDIR}/cmake/ttkBase"
  )

configure_file(core/base/TTKBaseConfig.cmake.in TTKBaseConfig.cmake @ONLY)
install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/TTKBaseConfig.cmake"
  DESTINATION
    "${CMAKE_INSTALL_LIBDIR}/cmake/ttkBase"
  )


# Paraview plugin
# ---------------

set(VTKWRAPPER_DIR "${CMAKE_CURRENT_LIST_DIR}/core/vtk/")

# VTK Wrappers
# ------------

if(TTK_BUILD_VTK_WRAPPERS)
  include(CMake/VTKModule.cmake)
  add_library(ttkAll INTERFACE)

  add_subdirectory(core/vtk)

  target_link_libraries(ttkAll INTERFACE ${TTK_ENABLED_MODULES})
  add_dependencies(ttkAll baseAll)

  install(
    TARGETS
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
      ${CMAKE_INSTALL_LIBDIR}/cmake/ttkVTK
    )

  configure_file(core/vtk/TTKVTKConfig.cmake.in TTKVTKConfig.cmake @ONLY)
  install(
    FILES
      "${CMAKE_CURRENT_BINARY_DIR}/TTKVTKConfig.cmake"
    DESTINATION
      ${CMAKE_INSTALL_LIBDIR}/cmake/ttkVTK
    )
endif()

if(TTK_BUILD_PARAVIEW_PLUGINS)
  include(CMake/ParaViewFilter.cmake)
  add_compile_definitions("TTK_BUILD_PARAVIEW_PLUGINS")

  # Install location
  if(NOT "$ENV{PV_PLUGIN_PATH}" STREQUAL "")
    set(TTK_INSTALL_PLUGIN_DIR
        $ENV{PV_PLUGIN_PATH}
        CACHE
        PATH
        "Directory where the ParaView plugin will be installed"
        )
  else()
    set(TTK_INSTALL_PLUGIN_DIR
        "${CMAKE_INSTALL_BINDIR}/plugins"
        CACHE
        PATH
        "Directory where the ParaView plugin will be installed"
        )
  endif()

  # set rpath for the plugin
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH};${TTK_INSTALL_PLUGIN_DIR}/${TTK_PLUGIN_SUBDIR}")

  add_subdirectory(paraview)
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
