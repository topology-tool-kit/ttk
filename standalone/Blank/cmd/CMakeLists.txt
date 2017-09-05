cmake_minimum_required(VERSION 3.1)

project(blankCmd)

add_executable(blankCmd main.cpp)
target_link_libraries(blankCmd ${VTK_LIBRARIES}
	ttkBlank ttkProgramBase)

