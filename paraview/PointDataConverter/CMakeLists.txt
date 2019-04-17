ttk_add_paraview_plugin(ttkPointDataConverter
  SOURCES
    ${VTKWRAPPER_DIR}/ttkPointDataConverter/ttkPointDataConverter.cpp
  PLUGIN_XML
    PointDataConverter.xml
  LINK
    common
    )

