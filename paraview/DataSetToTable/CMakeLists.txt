ttk_add_paraview_plugin(ttkDataSetToTable
  SOURCES
    ${VTKWRAPPER_DIR}/ttkDataSetToTable/ttkDataSetToTable.cpp
  PLUGIN_XML
    DataSetToTable.xml
  LINK
    common
    )

