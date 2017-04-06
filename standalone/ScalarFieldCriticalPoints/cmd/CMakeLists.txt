cmake_minimum_required(VERSION 2.4)

# name of the project (dummy)
project(scalarFieldCriticalPointsCmd)

# command-line build-mode
option(withVTKCMD "Enable VTK CMD support" true)
set(withVTKCMD true)

# path to the common code base
find_path(TTK_DIR ttk.cmake
  PATHS
    core/
    ../core/
    ../../core/
    ../../../core/
)
include(${TTK_DIR}/ttk.cmake)

# include the different packages required for the application
ttk_add_vtkWrapper_package(vtkScalarFieldCriticalPoints)
ttk_add_vtkWrapper_package(vtkProgramBase)

# add the local cpp files 
ttk_add_source(main.cpp)

# wrap-up the build
ttk_wrapup_binary(scalarFieldCriticalPointsCmd)
