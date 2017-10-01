ttk_add_paraview_plugin(ttkLDistance
	SOURCES ${VTKWRAPPER_DIR}/ttkLDistance/ttkLDistance.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML LDistance.xml
	LINK lDistance ttkTriangulation)

