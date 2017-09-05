add_paraview_plugin(ttkPVContourForests 1.0
	SERVER_MANAGER_XML ContourForests.xml
	SERVER_MANAGER_SOURCES ttkPVContourForests.cpp)
	#SERVER_MANAGER_SOURCES ../../core/vtkWrappers/ttkContourForests/ttkContourForests.cpp)
target_link_libraries(ttkPVContourForests PUBLIC ttkContourForests)
install(TARGETS ttkPVContourForests DESTINATION plugin)

