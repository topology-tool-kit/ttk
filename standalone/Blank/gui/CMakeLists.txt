cmake_minimum_required(VERSION 3.1)

project(blankGui)

find_package(VTK REQUIRED)
include_directories(${VTK_INCLUDE_DIRS})

add_executable(blankGui main.cpp)
target_link_libraries(blankGui ${VTK_LIBRARIES}
	ttkBlank ttkProgramBase ttkUserInterfaceBase)

