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
  )

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
  # Keep TTK libraries in the TopologyToolKit folder
  # foreach(TTK_MODULE ${TTK_MODULES})
  #   install(
  #     TARGETS
  #       ${TTK_MODULE}
  #     DESTINATION
  #       "${TTK_INSTALL_PLUGIN_DIR}/TopologyToolKit"
  #     )
  # endforeach()
  # set_target_properties(${TTK_PLUGINS}
  #   PROPERTIES INSTALL_RPATH
  #     "${TTK_INSTALL_PLUGIN_DIR}/TopologyToolKit"
  #   )

  # Install the TopologyToolKit library
  install(
    TARGETS
      TTKPVTargets
    DESTINATION
      "${TTK_INSTALL_PLUGIN_DIR}"
    )
endif()
