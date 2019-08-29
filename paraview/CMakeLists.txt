include(functions.cmake)

option(TTK_BUILD_STANDALONE_PARAVIEW_PLUGINS "Build one plugin per TTK filter" OFF)
mark_as_advanced(TTK_BUILD_STANDALONE_PARAVIEW_PLUGINS)

find_package(ParaView REQUIRED)
include(${PARAVIEW_USE_FILE})
include_directories(${PARAVIEW_INCLUDE_DIRS})

if (MSVC)
  set(TTK_INSTALL_PLUGIN_DIR
    "$ENV{PARAVIEW_BIN}/plugins"
    CACHE PATH
    "Directory where TTK ParaView plugins will be installed")
else()
  set(TTK_INSTALL_PLUGIN_DIR
    "${VTK_INSTALL_PREFIX}/lib/plugins"
    CACHE PATH
    "Directory where TTK ParaView plugins will be installed")
endif()

set(VTKWRAPPER_DIR "${CMAKE_CURRENT_LIST_DIR}/../core/vtk/")

# Build the ttkTriangulation library against the ParaView VTK
add_library(ttkPVTriangulation
  STATIC
    ${VTKWRAPPER_DIR}/ttkTriangulation/ttkTriangulation.cpp
    ${VTKWRAPPER_DIR}/ttkTriangulation/ttkWrapper.cpp
    )
target_link_libraries(ttkPVTriangulation PUBLIC ${VTK_LIBRARIES} triangulation)
set_target_properties(ttkPVTriangulation PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
target_include_directories(ttkPVTriangulation
  PUBLIC
    $<BUILD_INTERFACE:${VTKWRAPPER_DIR}/ttkTriangulation>
    )
if (MSVC)
  target_compile_definitions(ttkPVTriangulation PUBLIC TTK_PLUGIN)
endif()

set(TTK_PV_XML     "" CACHE INTERNAL "")
set(TTK_PV_SOURCES "" CACHE INTERNAL "")
set(TTK_PV_LINKS   "" CACHE INTERNAL "")

# traverse sub-directories
file(GLOB PARAVIEW_PLUGIN_DIRS *)
foreach(PARAVIEW_PLUGIN ${PARAVIEW_PLUGIN_DIRS})
  if (IS_DIRECTORY ${PARAVIEW_PLUGIN})
    if (NOT ${PARAVIEW_PLUGIN} MATCHES "patch" )
      add_subdirectory(${PARAVIEW_PLUGIN})
    endif()
  endif()
endforeach()

if (NOT TTK_BUILD_STANDALONE_PARAVIEW_PLUGINS)

  set_property(GLOBAL PROPERTY TTK_VTKWRAPPER_DIRS    "")
  # Build the ParaView plugin, if we're building with ParaView
  # and this wrapper library has a plugin
  add_paraview_plugin(TopologyToolKit ${PROJECT_VERSION}
    SERVER_MANAGER_XML
      ${TTK_PV_XML}
    SERVER_MANAGER_SOURCES
      ${TTK_PV_SOURCES}
      ${TTK_TRIANGULATION_SRCS}
      )

  if (MSVC)
    target_compile_definitions(TopologyToolKit PUBLIC TTK_PLUGIN)
  endif()

  # dependency on wrappers
  file(GLOB VTKWRAPPER_DIRS ${VTKWRAPPER_DIR}/*)
  foreach(VTKWRAPPER ${VTKWRAPPER_DIRS})
    if (IS_DIRECTORY ${VTKWRAPPER})
      list(APPEND TTK_VTKWRAPPER_DIRS ${VTKWRAPPER})
    endif()
  endforeach()

  target_include_directories(TopologyToolKit
    PUBLIC
      ${TTK_VTKWRAPPER_DIRS}
      )

  target_link_libraries(TopologyToolKit
    PUBLIC
      ${VTK_LIBRARIES}
      ${TTK_PV_LINKS}
      ttkPVTriangulation
      )

  install(TARGETS TopologyToolKit DESTINATION ${TTK_INSTALL_PLUGIN_DIR})
endif()
