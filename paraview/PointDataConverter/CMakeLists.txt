ttk_add_paraview_plugin(ttkPointDataConverter
	SOURCES ${VTKWRAPPER_DIR}/ttkPointDataConverter/ttkPointDataConverter.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML PointDataConverter.xml
	LINK common)

