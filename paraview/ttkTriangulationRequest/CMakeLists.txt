ttk_add_paraview_plugin(ttkTriangulationRequest
	SOURCES ${VTKWRAPPER_DIR}/ttkTriangulationRequest/ttkTriangulationRequest.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML TriangulationRequest.xml
	LINK ttkTriangulation)

