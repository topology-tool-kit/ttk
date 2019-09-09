ttk_add_paraview_plugin(ttkCellDataConverter
  SOURCES
    ${VTKWRAPPER_DIR}/ttkCellDataConverter/ttkCellDataConverter.cpp
  PLUGIN_XML
    CellDataConverter.xml
  LINK
    common
    )

