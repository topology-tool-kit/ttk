set(VTKWRAPPER_DIR "${CMAKE_CURRENT_LIST_DIR}/../core/vtk/")

# Fill TTK_MODULES, TTK_MODULES_FILES
# These lists are then used by singlePlugin/CMakeLists.txt
# The goal is to create a single plugin for parview
set(TTK_XMLS             "" CACHE INTERNAL "")
set(TTK_MODULES          "" CACHE INTERNAL "")
set(TTK_VTK_MODULE_FILES "" CACHE INTERNAL "")
file(GLOB PARAVIEW_PLUGIN_DIRS *)
foreach(PARAVIEW_PLUGIN ${PARAVIEW_PLUGIN_DIRS})
  if(IS_DIRECTORY ${PARAVIEW_PLUGIN})
    if(NOT ${PARAVIEW_PLUGIN} MATCHES "patch"
        AND NOT ${PARAVIEW_PLUGIN} MATCHES "singlePlugin")
      include(${PARAVIEW_PLUGIN}/TTKFilter.cmake)
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
    TTK_PLUGINS
  REQUIRES_MODULES
    TTK_REQUIRED
  )

foreach(REQ IN LISTS TTK_REQUIRED)
  if(NOT TARGET "${REQ}")
    message(FATAL_ERROR "Missing required module: ${REQ}")
  endif()
endforeach()

paraview_plugin_build(
  TARGET
    TTKPVTargets
  PLUGINS
    ${TTK_PLUGINS}
  AUTOLOAD
    ${TTK_PLUGINS}
  RUNTIME_DESTINATION
    "${CMAKE_INSTALL_BINDIR}"
  LIBRARY_DESTINATION
    "${CMAKE_INSTALL_LIBDIR}"
  LIBRARY_SUBDIRECTORY
    "${PARAVIEW_PLUGIN_SUBDIR}"
  )

if(NOT "${TTK_INSTALL_PLUGIN_DIR}" STREQUAL "")
  install(
    TARGETS
      TTKPVTargets
    DESTINATION
      "${TTK_INSTALL_PLUGIN_DIR}"
    )
endif()

