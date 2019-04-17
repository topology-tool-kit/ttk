ttk_add_paraview_plugin(ttkCellDataSelector
  SOURCES
    ${VTKWRAPPER_DIR}/ttkCellDataSelector/ttkCellDataSelector.cpp
  PLUGIN_XML
    CellDataSelector.xml
  LINK
    common
    )

