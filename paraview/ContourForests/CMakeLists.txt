ttk_add_paraview_plugin(ttkContourForests
  SOURCES
    ${VTKWRAPPER_DIR}/ttkContourForests/ttkContourForests.cpp
  PLUGIN_XML
    ContourForests.xml
  LINK
    geometry
    contourForests
    contourForestsTree
    )

