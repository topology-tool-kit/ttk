ttk_add_vtk_library(ttkBottleneckDistance
	SOURCES ttkBottleneckDistance.cpp
	HEADERS ttkBottleneckDistance.h
	PLUGIN_VERSION 1.0
	PLUGIN_XML BottleneckDistance.xml
	LINK bottleneckDistance ttkTriangulation)

