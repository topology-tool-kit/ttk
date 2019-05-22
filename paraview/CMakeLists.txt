include(functions.cmake)

if(MSVC)
  set(TTK_INSTALL_PLUGIN_DIR
    "$ENV{PARAVIEW_BIN}/plugins"
    CACHE
    PATH
    "Directory where TTK ParaView plugins will be installed")
else()
  set(TTK_INSTALL_PLUGIN_DIR
    "${VTK_INSTALL_PREFIX}/lib/plugins"
    CACHE
    PATH
    "Directory where TTK ParaView plugins will be installed")
endif()

set(VTKWRAPPER_DIR "${CMAKE_CURRENT_LIST_DIR}/../core/vtk/")

include(GNUInstallDirs)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")

# Fill TTK_MODULES, TTK_MODULES_FILES
# These lists are then used by singlePlugin/CMakeLists.txt
set(TTK_MODULES      "" CACHE INTERNAL "")
set(TTK_MODULE_FILES "" CACHE INTERNAL "")
file(GLOB PARAVIEW_PLUGIN_DIRS *)
foreach(PARAVIEW_PLUGIN ${PARAVIEW_PLUGIN_DIRS})
  if(IS_DIRECTORY ${PARAVIEW_PLUGIN})
    if(NOT ${PARAVIEW_PLUGIN} MATCHES "patch"
        AND NOT ${PARAVIEW_PLUGIN} MATCHES "singlePlugin")
      include(${PARAVIEW_PLUGIN}/TTKfilter.cmake)
    endif()
  endif()
endforeach()

set("PARAVIEW_PLUGIN_ENABLE_TopologyToolKit"
  ON
  CACHE
  BOOL
  "Build TTK TopologyToolKit filter"
  FORCE
  )

paraview_plugin_scan(
  PLUGIN_FILES
  "singlePlugin/paraview.plugin"
  PROVIDES_PLUGINS
    plugins
  )

paraview_plugin_build(
  RUNTIME_DESTINATION
    "${CMAKE_INSTALL_BINDIR}"
  LIBRARY_DESTINATION
    "${CMAKE_INSTALL_LIBDIR}"
  LIBRARY_SUBDIRECTORY
    ""
  PLUGINS
    ${plugins}
    )
