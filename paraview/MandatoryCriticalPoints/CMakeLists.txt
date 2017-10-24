ttk_add_paraview_plugin(ttkMandatoryCriticalPoints
	SOURCES ${VTKWRAPPER_DIR}/ttkMandatoryCriticalPoints/ttkMandatoryCriticalPoints.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML MandatoryCriticalPoints.xml
	LINK mandatoryCriticalPoints)

