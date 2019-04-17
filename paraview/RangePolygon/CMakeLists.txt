ttk_add_paraview_plugin(ttkRangePolygon
  SOURCES
    ${VTKWRAPPER_DIR}/ttkRangePolygon/ttkRangePolygon.cpp
  PLUGIN_XML
    RangePolygon.xml
  LINK
    scalarFieldSmoother
    )

