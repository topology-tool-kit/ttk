ttk_add_paraview_plugin(ttkBlank
	SOURCES ${VTKWRAPPER_DIR}/ttkBlank/ttkBlank.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML Blank.xml
	LINK blank ttkTriangulation)

