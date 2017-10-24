ttk_add_paraview_plugin(ttkCellDataConverter
	SOURCES ${VTKWRAPPER_DIR}/ttkCellDataConverter/ttkCellDataConverter.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML CellDataConverter.xml
	LINK common)

