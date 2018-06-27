ttk_add_paraview_plugin(ttkTopologicalCompressionReader
	SOURCES ${VTKWRAPPER_DIR}/ttkTopologicalCompressionReader/ttkTopologicalCompressionReader.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML TopologicalCompressionReader.xml
  	LINK topologicalCompression)
