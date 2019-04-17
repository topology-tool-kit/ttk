ttk_add_paraview_plugin(ttkIntegralLines
  SOURCES
    ${VTKWRAPPER_DIR}/ttkIntegralLines/ttkIntegralLines.cpp
  PLUGIN_XML
    IntegralLines.xml
  LINK
    integralLines
    )

