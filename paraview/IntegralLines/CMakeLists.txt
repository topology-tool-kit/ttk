ttk_add_paraview_plugin(ttkIntegralLines
	SOURCES ${VTKWRAPPER_DIR}/ttkIntegralLines/ttkIntegralLines.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML IntegralLines.xml
	LINK integralLines)

