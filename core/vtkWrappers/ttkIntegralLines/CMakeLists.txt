ttk_add_vtk_library(ttkIntegralLines
	SOURCES ttkIntegralLines.cpp
	HEADERS ttkIntegralLines.h
	PLUGIN_VERSION 1.0
	PLUGIN_XML IntegralLines.xml
	LINK integralLines ttkTriangulation)

