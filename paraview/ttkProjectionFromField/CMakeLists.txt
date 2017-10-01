ttk_add_paraview_plugin(ttkProjectionFromField
	SOURCES ${VTKWRAPPER_DIR}/ttkProjectionFromField/ttkProjectionFromField.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML ProjectionFromField.xml
	LINK common)

