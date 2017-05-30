cmake_minimum_required(VERSION 2.4)

project(IdentifierRandomizer)

set(TTK_BUILD_MODE "ParaView")

# path to the common code base
find_path(TTK_DIR ttk.cmake
  PATHS
    core/
    ../core/
    ../../core/
    ../../../core/
    sandbox/core/
    sandbox/../core/
    sandbox/../../core/
    sandbox/../../../core/
)
include(${TTK_DIR}/ttk.cmake)

# include the vtk wrapper for your TTK package 
ttk_add_vtkWrapper_package(ttkIdentifierRandomizer)

# Add the plugin
ttk_wrapup_paraview_plugin(
  IdentifierRandomizer
  1.0)
