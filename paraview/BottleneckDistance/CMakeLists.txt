cmake_minimum_required(VERSION 2.4)

project(BottleneckDistance)

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
ttk_add_vtkWrapper_package(ttkBottleneckDistance)

# Add the plugin
ttk_wrapup_paraview_plugin(
  BottleneckDistance
  1.0)

