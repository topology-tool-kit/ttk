ttk_add_paraview_plugin(ttkMeshSubdivision
	SOURCES ${VTKWRAPPER_DIR}/ttkMeshSubdivision/ttkMeshSubdivision.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML MeshSubdivision.xml
	LINK common)

