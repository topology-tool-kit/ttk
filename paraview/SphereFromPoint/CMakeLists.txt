ttk_add_paraview_plugin(ttkSphereFromPoint
	SOURCES ${VTKWRAPPER_DIR}/ttkSphereFromPoint/ttkSphereFromPoint.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML SphereFromPoint.xml
	LINK common)

