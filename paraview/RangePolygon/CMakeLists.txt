ttk_add_paraview_plugin(ttkRangePolygon
	SOURCES ${VTKWRAPPER_DIR}/ttkRangePolygon/ttkRangePolygon.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML RangePolygon.xml
	LINK scalarFieldSmoother)

