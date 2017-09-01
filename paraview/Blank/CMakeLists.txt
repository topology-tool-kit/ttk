cmake_minimum_required(VERSION 3.1)
project(Blank)

add_paraview_plugin(Blank 1.0
	SERVER_MANAGER_XML Blank.xml
	SERVER_MANAGER_SOURCES Blank.cpp)

target_link_libraries(Blank ttk::vtk::Blank)

install(TARGETS Blank DESTINATION bin)

