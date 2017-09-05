ttk_add_vtk_library(ttkContourForests ttkContourForests.cpp
	HEADERS ttkContourForests.h)
target_link_libraries(ttkContourForests PUBLIC geometry contourForests
	contourForestsTree ttkTriangulation)

if (TTK_BUILD_PARAVIEW_PLUGINS)
	add_paraview_plugin(ttkContourForestsPlugin 1.0
		SERVER_MANAGER_XML ContourForests.xml
		SERVER_MANAGER_SOURCES ttkContourForests.cpp)
	target_link_libraries(ttkContourForestsPlugin PUBLIC geometry contourForests
		contourForestsTree ttkTriangulation)
	install(TARGETS ttkContourForestsPlugin DESTINATION plugin)
endif()

