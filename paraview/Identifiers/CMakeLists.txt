ttk_add_paraview_plugin(ttkIdentifiers
  SOURCES
    ${VTKWRAPPER_DIR}/ttkIdentifiers/ttkIdentifiers.cpp
  PLUGIN_XML
    Identifiers.xml
  LINK
    common
    )

