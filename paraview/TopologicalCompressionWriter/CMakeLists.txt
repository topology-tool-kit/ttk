ttk_add_paraview_plugin(ttkTopologicalCompressionWriter
	SOURCES ${VTKWRAPPER_DIR}/ttkTopologicalCompressionWriter/ttkTopologicalCompressionWriter.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML TopologicalCompressionWriter.xml
  	LINK topologicalCompression)
