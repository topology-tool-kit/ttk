ttk_add_paraview_plugin(
  ttkSurfaceQuadrangulation
  SOURCES
  ${VTKWRAPPER_DIR}/ttkSurfaceQuadrangulation/ttkSurfaceQuadrangulation.cpp
  PLUGIN_XML
  SurfaceQuadrangulation.xml
  LINK
  surfaceQuadrangulation)
