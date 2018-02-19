cmake_minimum_required(VERSION 3.2)

project(manifoldCheckCmd)

set(CMAKE_SKIP_BUILD_RPATH TRUE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE) 
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/ttk/")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_executable(manifoldCheckCmd main.cpp)
target_link_libraries(manifoldCheckCmd ${VTK_LIBRARIES}
	ttkManifoldCheck ttkProgramBase)

install(TARGETS manifoldCheckCmd RUNTIME DESTINATION ${TTK_INSTALL_BINARY_DIR})
