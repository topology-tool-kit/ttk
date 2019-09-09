ttk_add_paraview_plugin(ttkProjectionFromField
  SOURCES
    ${VTKWRAPPER_DIR}/ttkProjectionFromField/ttkProjectionFromField.cpp
  PLUGIN_XML
    ProjectionFromField.xml
  LINK
    common
    )

