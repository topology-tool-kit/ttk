ttk_add_paraview_plugin(ttkReebSpace
	SOURCES ${VTKWRAPPER_DIR}/ttkReebSpace/ttkReebSpace.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML ReebSpace.xml
	LINK reebSpace ttkTriangulation)

