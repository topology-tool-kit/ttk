include(functions.cmake)

add_subdirectory(base)

if(TTK_BUILD_PARAVIEW_PLUGINS)
  # If paraview plugin is on, modules will be automatically built by parview
  include(vtk/functions.cmake)
  include(vtk/ttkTriangulation/TTKWrapper.cmake)
elseif(TTK_BUILD_VTK_WRAPPERS)
  add_subdirectory(vtk)
endif()
