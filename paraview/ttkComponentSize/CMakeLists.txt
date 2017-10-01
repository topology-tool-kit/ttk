ttk_add_paraview_plugin(ttkComponentSize
	SOURCES ${VTKWRAPPER_DIR}/ttkComponentSize/ttkComponentSize.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML ComponentSize.xml
	LINK common)

