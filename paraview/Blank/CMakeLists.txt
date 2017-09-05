add_paraview_plugin(ttkPVBlank 1.0
	SERVER_MANAGER_XML Blank.xml)
target_link_libraries(ttkPVBlank PUBLIC ttk::vtk::Blank)
install(TARGETS ttkPVBlank DESTINATION plugin)

