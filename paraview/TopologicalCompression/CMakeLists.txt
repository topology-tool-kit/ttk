ttk_add_paraview_plugin(ttkTopologicalCompression
  SOURCES
    ${VTKWRAPPER_DIR}/ttkTopologicalCompression/ttkTopologicalCompression.cpp
  PLUGIN_XML
    TopologicalCompression.xml
  LINK
    topologicalCompression
    )
