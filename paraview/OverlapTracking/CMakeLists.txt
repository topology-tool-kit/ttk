ttk_add_paraview_plugin(ttkOverlapTracking
	SOURCES ${VTKWRAPPER_DIR}/ttkOverlapTracking/ttkOverlapTracking.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML OverlapTracking.xml
	LINK overlapTracking)