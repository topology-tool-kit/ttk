cmake_minimum_required(VERSION 2.4)

project(LDistance)

set(TTK_BUILD_MODE "ParaView")

# path to the common code base
find_path(TTK_DIR ttk.cmake
  PATHS
    core/
    ../core/
    ../../core/
    ../../../core/
)
include(${TTK_DIR}/ttk.cmake)

# include the vtk wrapper for your TTK package 
ttk_add_vtkWrapper_package(ttkLDistance)

# Add the plugin
ttk_wrapup_paraview_plugin(
  LDistance
  1.0)

