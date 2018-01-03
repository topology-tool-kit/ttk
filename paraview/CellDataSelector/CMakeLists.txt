ttk_add_paraview_plugin(ttkCellDataSelector
	SOURCES ${VTKWRAPPER_DIR}/ttkCellDataSelector/ttkCellDataSelector.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML CellDataSelector.xml
	LINK common)

