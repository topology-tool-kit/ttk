ttk_add_paraview_plugin(ttkFTRGraph
  SOURCES
    ${VTKWRAPPER_DIR}/ttkFTRGraph/ttkFTRGraph.cpp
  PLUGIN_XML
    FTRGraph.xml
  LINK
    ftrGraph
    )

