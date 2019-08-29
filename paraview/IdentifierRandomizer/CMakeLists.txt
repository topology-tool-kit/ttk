ttk_add_paraview_plugin(ttkIdentifierRandomizer
  SOURCES
    ${VTKWRAPPER_DIR}/ttkIdentifierRandomizer/ttkIdentifierRandomizer.cpp
  PLUGIN_XML
    IdentifierRandomizer.xml
  LINK
    common
    )

