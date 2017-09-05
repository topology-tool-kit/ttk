ttk_add_vtk_library(ttkFiber
	SOURCES ttkFiber.cpp
	HEADERS ttkFiber.h
	PLUGIN_VERSION 1.0
	PLUGIN_XML Fiber.xml
	LINK common ttkTriangulation)

