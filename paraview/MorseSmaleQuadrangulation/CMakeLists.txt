ttk_add_paraview_plugin(
  ttkMorseSmaleQuadrangulation
  SOURCES
  ${VTKWRAPPER_DIR}/ttkMorseSmaleQuadrangulation/ttkMorseSmaleQuadrangulation.cpp
  PLUGIN_XML
  MorseSmaleQuadrangulation.xml
  LINK
  morseSmaleQuadrangulation)
