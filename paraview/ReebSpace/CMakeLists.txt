ttk_add_paraview_plugin(ttkReebSpace
  SOURCES
    ${VTKWRAPPER_DIR}/ttkReebSpace/ttkReebSpace.cpp
  PLUGIN_XML
    ReebSpace.xml
  LINK
    reebSpace
    )

