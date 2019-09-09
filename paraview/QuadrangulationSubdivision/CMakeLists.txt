ttk_add_paraview_plugin(
  ttkQuadrangulationSubdivision
  SOURCES
  ${VTKWRAPPER_DIR}/ttkQuadrangulationSubdivision/ttkQuadrangulationSubdivision.cpp
  PLUGIN_XML
  QuadrangulationSubdivision.xml
  LINK
  quadrangulationSubdivision)
