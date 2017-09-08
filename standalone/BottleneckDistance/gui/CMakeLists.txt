cmake_minimum_required(VERSION 3.5)

project(bottleneckDistanceGui)

add_executable(bottleneckDistanceGui main.cpp)
target_link_libraries(bottleneckDistanceGui  ${VTK_LIBRARIES}
	ttkBottleneckDistance ttkUserInterfaceBase)
install(TARGETS bottleneckDistanceGui RUNTIME DESTINATION bin/ttk)

