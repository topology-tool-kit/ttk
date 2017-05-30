cmake_minimum_required(VERSION 2.4)

project(ScalarFieldCriticalPoints)

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

# locate paraview
find_package(ParaView REQUIRED)
include(${PARAVIEW_USE_FILE})

# include the vtk wrapper for the scalar field smoothing package
# HERE: change to whatever wrapper you want to use
ttk_add_vtkWrapper_package(ttkScalarFieldCriticalPoints)

# Add the plugin
ttk_wrapup_paraview_plugin(
  ScalarFieldCriticalPoints
  1.0)
