ttk_add_paraview_plugin(ttkManifoldCheck
	SOURCES ${VTKWRAPPER_DIR}/ttkManifoldCheck/ttkManifoldCheck.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML ManifoldCheck.xml
	LINK manifoldCheck)

