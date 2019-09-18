set(VTKWRAPPER_DIR "${CMAKE_CURRENT_LIST_DIR}/../core/vtk/")

# Fill TTK_MODULES, TTK_MODULES_FILES
# These lists are then used by singlePlugin/CMakeLists.txt
set(TTK_XMLS              "" CACHE INTERNAL "")
set(TTK_MODULES           "" CACHE INTERNAL "")
set(TTK_VTK_MODULE_FILESS "" CACHE INTERNAL "")
file(GLOB PARAVIEW_PLUGIN_DIRS *)
foreach(PARAVIEW_PLUGIN ${PARAVIEW_PLUGIN_DIRS})
  if(IS_DIRECTORY ${PARAVIEW_PLUGIN})
    if(NOT ${PARAVIEW_PLUGIN} MATCHES "patch"
        AND NOT ${PARAVIEW_PLUGIN} MATCHES "singlePlugin")
      include(${PARAVIEW_PLUGIN}/TTKFilter.cmake)
    endif()
  endif()
endforeach()

if(NOT "$ENV{PV_PLUGIN_PATH}" STREQUAL "")
  set(TTK_INSTALL_PLUGIN_DIR
    $ENV{PV_PLUGIN_PATH}
    CACHE
    PATH
      "Directory where the ParaView plugin will be installed"
    )
else()
  message(WARNING "Please, set the PV_PLUGIN_PATH environement variable, the TTK ParaView plugin will be built but not installed.")
endif()

include(GNUInstallDirs)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")

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

if(DEFINED TTK_INSTALL_PLUGIN_DIR)
  # each filter library
  foreach(TTK_MODULE ${TTK_MODULES})
    install(
      TARGETS
        ${TTK_MODULE}
      DESTINATION
        "${TTK_INSTALL_PLUGIN_DIR}/TopologyToolKit"
      )
  endforeach()
  # the main plugin
  set_target_properties(${TTK_PLUGINS}
    PROPERTIES INSTALL_RPATH
      "${TTK_INSTALL_PLUGIN_DIR}/TopologyToolKit"
    )
  install(
    TARGETS
      ${TTK_PLUGINS}
    DESTINATION
      "${TTK_INSTALL_PLUGIN_DIR}"
    )
endif()
