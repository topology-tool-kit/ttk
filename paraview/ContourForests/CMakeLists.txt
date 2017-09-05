add_paraview_plugin(ttkPVContourForests 1.0
	SERVER_MANAGER_XML ContourForests.xml)
target_link_libraries(ttkPVContourForests PUBLIC ttk::vtk::ContourForests)
install(TARGETS ttkPVContourForests DESTINATION plugin)

