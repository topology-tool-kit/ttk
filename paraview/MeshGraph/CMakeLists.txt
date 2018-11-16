ttk_add_paraview_plugin(ttkMeshGraph
	SOURCES ${VTKWRAPPER_DIR}/ttkMeshGraph/ttkMeshGraph.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML MeshGraph.xml
	LINK meshGraph)

