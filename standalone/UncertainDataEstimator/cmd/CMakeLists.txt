cmake_minimum_required(VERSION 3.2)

project(uncertainDataEstimatorCmd)

add_executable(uncertainDataEstimatorCmd main.cpp Editor.cpp)
target_link_libraries(uncertainDataEstimatorCmd ttkUncertainDataEstimator ttkProgramBase)
target_include_directories(uncertainDataEstimatorCmd PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
install(TARGETS uncertainDataEstimatorCmd RUNTIME DESTINATION bin/ttk)

