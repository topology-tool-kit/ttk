ttk_add_paraview_plugin(ttkPointDataSelector
	SOURCES ${VTKWRAPPER_DIR}/ttkPointDataSelector/ttkPointDataSelector.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML PointDataSelector.xml
	LINK common)

