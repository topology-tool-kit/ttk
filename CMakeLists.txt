cmake_minimum_required(VERSION 3.1)
project(ttk)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_subdirectory(core)
add_subdirectory(standalone)
add_subdirectory(paraview)

