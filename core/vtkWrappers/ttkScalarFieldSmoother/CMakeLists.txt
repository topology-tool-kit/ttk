ttk_add_vtk_library(ttkScalarFieldSmoother
	SOURCES ttkScalarFieldSmoother.cpp
	HEADERS ttkScalarFieldSmoother.h
	PLUGIN_VERSION 1.0
	PLUGIN_XML ScalarFieldSmoother.xml
	LINK scalarFieldSmoother ttkTriangulation)

