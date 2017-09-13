cmake_minimum_required(VERSION 3.2)

project(continuousScatterPlotCmd)

add_executable(continuousScatterPlotCmd main.cpp)
target_link_libraries(continuousScatterPlotCmd ttkContinuousScatterPlot
	ttkProgramBase)

install(TARGETS continuousScatterPlotCmd RUNTIME DESTINATION ${TTK_INSTALL_BINARY_DIR})

