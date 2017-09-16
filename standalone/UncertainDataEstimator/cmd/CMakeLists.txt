cmake_minimum_required(VERSION 3.2)

project(uncertainDataEstimatorCmd)

set(CMAKE_SKIP_BUILD_RPATH TRUE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE) 
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/ttk/")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_executable(uncertainDataEstimatorCmd main.cpp Editor.cpp)
target_link_libraries(uncertainDataEstimatorCmd ttkUncertainDataEstimator ttkProgramBase)
target_include_directories(uncertainDataEstimatorCmd PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
install(TARGETS uncertainDataEstimatorCmd RUNTIME DESTINATION ${TTK_INSTALL_BINARY_DIR})

