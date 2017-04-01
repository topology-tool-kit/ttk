cmake_minimum_required(VERSION 2.4)

# name of the project (persistenceDiagramGui)
project(persistenceDiagramGui)

# GUI build-mode
option(withVTKGUI "Enable VTK GUI support" true)
set(withVTKGUI true)

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
ttk_add_vtkWrapper_package(vtkPersistenceDiagram)

# add the local cpp files
ttk_add_source(main.cpp)

ttk_wrapup_binary(persistenceDiagramGui)
