ttk_add_paraview_plugin(ttkDistanceField
	SOURCES ${VTKWRAPPER_DIR}/ttkDistanceField/ttkDistanceField.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML DistanceField.xml
	LINK distanceField ttkTriangulation)

