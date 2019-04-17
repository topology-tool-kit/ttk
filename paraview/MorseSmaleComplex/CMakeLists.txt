ttk_add_paraview_plugin(ttkMorseSmaleComplex
  SOURCES
    ${VTKWRAPPER_DIR}/ttkMorseSmaleComplex/ttkMorseSmaleComplex.cpp
  PLUGIN_XML
    MorseSmaleComplex.xml
  LINK
    morseSmaleComplex
    )

