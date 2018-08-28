ttk_add_paraview_plugin(ttkManifoldLearning
	SOURCES ${VTKWRAPPER_DIR}/ttkManifoldLearning/ttkManifoldLearning.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML ManifoldLearning.xml
	LINK manifoldLearning)

