cmake_minimum_required(VERSION 2.4)

# name of the project (dummy)
project(geometrySmootherCmd)

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
ttk_add_baseCode_package(morseSmaleComplex)
ttk_add_baseCode_package(persistenceCurve)
ttk_add_baseCode_package(persistenceDiagram)
ttk_add_baseCode_package(scalarFieldSmoother)
ttk_add_baseCode_package(topologicalSimplification)

# add the local cpp files 
ttk_add_source(main.cpp)

# wrap-up the build
ttk_wrapup_binary(geometrySmootherCmd)
