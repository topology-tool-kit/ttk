ttk_add_vtk_library(ttkRangePolygon
	SOURCES ttkRangePolygon.cpp
	HEADERS ttkRangePolygon.h
	PLUGIN_VERSION 1.0
	PLUGIN_XML RangePolygon.xml
	LINK scalarFieldSmoother ttkTriangulation)

