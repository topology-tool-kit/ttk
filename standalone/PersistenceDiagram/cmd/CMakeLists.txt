cmake_minimum_required(VERSION 3.5)

project(persistenceDiagramCmd)
add_executable(persistenceDiagramCmd main.cpp)
target_link_libraries(persistenceDiagramCmd ttkPersistenceDiagram ttkProgramBase)

