ttk_add_vtk_library(ttkMandatoryCriticalPoints
	SOURCES ttkMandatoryCriticalPoints.cpp
	HEADERS ttkMandatoryCriticalPoints.h
	PLUGIN_VERSION 1.0
	PLUGIN_XML MandatoryCriticalPoints.xml
	LINK mandatoryCriticalPoints ttkTriangulation)

