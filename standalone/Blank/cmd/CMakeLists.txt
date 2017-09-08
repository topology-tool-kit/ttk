cmake_minimum_required(VERSION 3.5)

project(blankCmd)

add_executable(blankCmd main.cpp)
target_link_libraries(blankCmd ${VTK_LIBRARIES}
	ttkBlank ttkProgramBase)
install(TARGETS blankCmd RUNTIME DESTINATION bin/ttk)

