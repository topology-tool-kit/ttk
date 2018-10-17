ttk_add_paraview_plugin(ttkImageMetrics
	SOURCES ${VTKWRAPPER_DIR}/ttkImageMetrics/ttkImageMetrics.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML ImageMetrics.xml
	LINK imageMetrics)

