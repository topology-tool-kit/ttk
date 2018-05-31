ttk_add_paraview_plugin(ttkCompare
	SOURCES ${VTKWRAPPER_DIR}/ttkCompare/ttkCompare.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML Compare.xml
	LINK compare)

