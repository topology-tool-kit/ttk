ttk_add_paraview_plugin(ttkBottleneckDistance
	SOURCES ${VTKWRAPPER_DIR}/ttkBottleneckDistance/ttkBottleneckDistance.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML BottleneckDistance.xml
	LINK bottleneckDistance)

