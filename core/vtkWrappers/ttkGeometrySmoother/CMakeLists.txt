ttk_add_vtk_library(ttkGeometrySmoother
	SOURCES ttkGeometrySmoother.cpp
	HEADERS ttkGeometrySmoother.h
	PLUGIN_VERSION 1.0
	PLUGIN_XML GeometrySmoother.xml
	LINK scalarFieldSmoother ttkTriangulation)

