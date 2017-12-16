ttk_add_paraview_plugin(ttkFieldSelector
	SOURCES ${VTKWRAPPER_DIR}/ttkFieldSelector/ttkFieldSelector.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML FieldSelector.xml
	LINK common)

