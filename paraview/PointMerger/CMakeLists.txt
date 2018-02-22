ttk_add_paraview_plugin(ttkPointMerger
	SOURCES ${VTKWRAPPER_DIR}/ttkPointMerger/ttkPointMerger.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML PointMerger.xml
	LINK pointMerger)

