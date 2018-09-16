ttk_add_paraview_plugin(ttkA3Refinement
	SOURCES ${VTKWRAPPER_DIR}/ttkA3Refinement/ttkA3Refinement.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML A3Refinement.xml
	LINK a3Refinement)

