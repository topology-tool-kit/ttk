ttk_add_paraview_plugin(ttkDataSetToTable
	SOURCES ${VTKWRAPPER_DIR}/ttkDataSetToTable/ttkDataSetToTable.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML DataSetToTable.xml
	LINK common)

