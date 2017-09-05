ttk_add_vtk_library(ttkFiberSurface
	SOURCES ttkFiberSurface.cpp
	HEADERS ttkFiberSurface.h
	PLUGIN_VERSION 1.0
	PLUGIN_XML FiberSurface.xml
	LINK fiberSurface ttkTriangulation)


