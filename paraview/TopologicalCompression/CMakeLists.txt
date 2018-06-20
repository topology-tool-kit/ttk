ttk_add_paraview_plugin(ttkTopologicalCompression
	SOURCES ${VTKWRAPPER_DIR}/ttkTopologicalCompression/ttkTopologicalCompression.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML TopologicalCompression.xml
	LINK topologicalCompression)
