ttk_add_paraview_plugin(ttkFiber
	SOURCES ${VTKWRAPPER_DIR}/ttkFiber/ttkFiber.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML Fiber.xml
	LINK common ttkTriangulation)

