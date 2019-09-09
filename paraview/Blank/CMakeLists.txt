ttk_add_paraview_plugin(ttkBlank
  SOURCES
    ${VTKWRAPPER_DIR}/ttkBlank/ttkBlank.cpp
  PLUGIN_XML
    Blank.xml
  LINK
    blank
    )

