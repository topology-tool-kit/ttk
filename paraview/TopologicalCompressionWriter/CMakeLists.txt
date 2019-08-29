ttk_add_paraview_plugin(ttkTopologicalCompressionWriter
  SOURCES
    ${VTKWRAPPER_DIR}/ttkTopologicalCompressionWriter/ttkTopologicalCompressionWriter.cpp
  PLUGIN_XML
    TopologicalCompressionWriter.xml
  LINK
    topologicalCompression
    )
