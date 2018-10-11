ttk_add_paraview_plugin(ttkTableDataSelector
	SOURCES ${VTKWRAPPER_DIR}/ttkTableDataSelector/ttkTableDataSelector.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML TableDataSelector.xml
	LINK common)

