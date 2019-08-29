ttk_add_paraview_plugin(ttkHarmonicField
  SOURCES
    ${VTKWRAPPER_DIR}/ttkHarmonicField/ttkHarmonicField.cpp
  PLUGIN_XML
    HarmonicField.xml
  LINK
    harmonicField
    )
