add_paraview_plugin(Blank 1.0
	SERVER_MANAGER_XML Blank.xml)

target_link_libraries(Blank PUBLIC ttk::vtk::Blank)

install(TARGETS Blank DESTINATION bin)

