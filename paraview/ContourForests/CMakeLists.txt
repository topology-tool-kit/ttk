add_paraview_plugin(ContourForests 1.0
	SERVER_MANAGER_XML ContourForests.xml)
target_link_libraries(ContourForests PUBLIC ttk::vtk::ContourForests)
install(TARGETS ContourForests DESTINATION plugin)

