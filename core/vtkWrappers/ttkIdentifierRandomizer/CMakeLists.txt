ttk_add_vtk_library(ttkIdentifierRandomizer
	SOURCES ttkIdentifierRandomizer.cpp
	HEADERS ttkIdentifierRandomizer.h
	PLUGIN_VERSION 1.0
	PLUGIN_XML IdentifierRandomizer.xml
	LINK common ttkTriangulation)

