cmake_minimum_required(VERSION 3.2)

project(ttkPersistenceDiagramsClusteringCmd)

set(CMAKE_SKIP_BUILD_RPATH TRUE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE) 
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/ttk/")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_executable(ttkPersistenceDiagramsClusteringCmd main.cpp )
target_link_libraries(ttkPersistenceDiagramsClusteringCmd ttkPersistenceDiagramsClustering ttkProgramBase)
target_include_directories(ttkPersistenceDiagramsClusteringCmd PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
install(TARGETS ttkPersistenceDiagramsClusteringCmd RUNTIME DESTINATION ${TTK_INSTALL_BINARY_DIR})

