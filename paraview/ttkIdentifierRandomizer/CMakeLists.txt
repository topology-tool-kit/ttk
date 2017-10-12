ttk_add_paraview_plugin(ttkIdentifierRandomizer
	SOURCES ${VTKWRAPPER_DIR}/ttkIdentifierRandomizer/ttkIdentifierRandomizer.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML IdentifierRandomizer.xml
	LINK common)

