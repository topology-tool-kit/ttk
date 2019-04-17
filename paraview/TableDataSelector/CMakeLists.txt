ttk_add_paraview_plugin(ttkTableDataSelector
  SOURCES
    ${VTKWRAPPER_DIR}/ttkTableDataSelector/ttkTableDataSelector.cpp
  PLUGIN_XML
    TableDataSelector.xml
  LINK
    common
    )

