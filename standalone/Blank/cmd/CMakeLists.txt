cmake_minimum_required(VERSION 3.1)

project(blankCmd)

find_package(VTK REQUIRED)
include_directories(${VTK_INCLUDE_DIRS})

add_executable(blankCmd main.cpp)
target_link_libraries(blankCmd ${VTK_LIBRARIES}
	ttk::vtk::Blank vtkProgramBase)

