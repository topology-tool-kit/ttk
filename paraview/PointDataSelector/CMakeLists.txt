ttk_add_paraview_plugin(ttkPointDataSelector
  SOURCES
    ${VTKWRAPPER_DIR}/ttkPointDataSelector/ttkPointDataSelector.cpp
  PLUGIN_XML
    PointDataSelector.xml
  LINK
    common
    )

