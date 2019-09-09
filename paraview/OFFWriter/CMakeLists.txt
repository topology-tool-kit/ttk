ttk_add_paraview_plugin(ttkOFFWriter
  SOURCES
    ${VTKWRAPPER_DIR}/ttkOFFWriter/ttkOFFWriter.cpp
  PLUGIN_XML
    OFFWriter.xml
    )
