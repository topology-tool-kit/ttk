ttk_add_vtk_library(ttkDistanceField
	SOURCES ttkDistanceField.cpp
	HEADERS ttkDistanceField.h
	PLUGIN_VERSION 1.0
	PLUGIN_XML DistanceField.xml
	LINK distanceField ttkTriangulation)

