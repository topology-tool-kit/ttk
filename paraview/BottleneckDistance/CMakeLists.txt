ttk_add_paraview_plugin(ttkBottleneckDistance
  SOURCES
    ${VTKWRAPPER_DIR}/ttkBottleneckDistance/ttkBottleneckDistance.cpp
  PLUGIN_XML
    BottleneckDistance.xml
  LINK
    bottleneckDistance
    )

