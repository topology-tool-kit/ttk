cmake_minimum_required(VERSION 3.5)

project(meshSubdivisionCmd)
add_executable(meshSubdivisionCmd main.cpp)
target_link_libraries(meshSubdivisionCmd ttkMeshSubdivision ttkProgramBase)

