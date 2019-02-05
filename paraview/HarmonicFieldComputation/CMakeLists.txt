ttk_add_paraview_plugin(
  ttkHarmonicFieldComputation
  SOURCES
  ${VTKWRAPPER_DIR}/ttkHarmonicFieldComputation/ttkHarmonicFieldComputation.cpp
  PLUGIN_XML
  HarmonicFieldComputation.xml
  LINK
  harmonicFieldComputation)
