ttk_add_paraview_plugin(ttkComponentSize
  SOURCES
    ${VTKWRAPPER_DIR}/ttkComponentSize/ttkComponentSize.cpp
  PLUGIN_XML
    ComponentSize.xml
  LINK
    common
    )

