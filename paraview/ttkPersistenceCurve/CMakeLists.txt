ttk_add_paraview_plugin(ttkPersistenceCurve
	SOURCES ${VTKWRAPPER_DIR}/ttkPersistenceCurve/ttkPersistenceCurve.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML PersistenceCurve.xml
	LINK persistenceCurve)

