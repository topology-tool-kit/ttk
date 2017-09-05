ttk_add_vtk_library(ttkContourForests
	SOURCES ttkContourForests.cpp
	HEADERS ttkContourForests.h
	PLUGIN_VERSION 1.0
	PLUGIN_XML ContourForests.xml
	LINK geometry contourForests contourForestsTree ttkTriangulation)

