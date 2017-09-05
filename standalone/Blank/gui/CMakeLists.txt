cmake_minimum_required(VERSION 3.1)

project(blankGui)

add_executable(blankGui main.cpp)
target_link_libraries(blankGui ${VTK_LIBRARIES}
	ttkBlank ttkProgramBase ttkUserInterfaceBase)

