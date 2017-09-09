cmake_minimum_required(VERSION 3.2)

project(bottleneckDistanceGui)

add_executable(bottleneckDistanceGui main.cpp)
target_link_libraries(bottleneckDistanceGui  ${VTK_LIBRARIES}
	ttkBottleneckDistance ttkUserInterfaceBase)
install(TARGETS bottleneckDistanceGui RUNTIME DESTINATION bin/ttk)

