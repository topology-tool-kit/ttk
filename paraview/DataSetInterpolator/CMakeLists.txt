ttk_add_paraview_plugin(ttkDataSetInterpolator
  SOURCES
    ${VTKWRAPPER_DIR}/ttkDataSetInterpolator/ttkDataSetInterpolator.cpp
  PLUGIN_XML
    DataSetInterpolator.xml
  LINK
    common
    )

