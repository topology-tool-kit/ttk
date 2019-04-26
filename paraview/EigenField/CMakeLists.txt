ttk_add_paraview_plugin(ttkEigenField
  SOURCES
    ${VTKWRAPPER_DIR}/ttkEigenField/ttkEigenField.cpp
  PLUGIN_XML
    EigenField.xml
  LINK
    eigenField
    )
