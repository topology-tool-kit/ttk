ttk_add_paraview_plugin(ttkDepthImageBasedGeometryApproximation
  SOURCES
    ${VTKWRAPPER_DIR}/ttkDepthImageBasedGeometryApproximation/ttkDepthImageBasedGeometryApproximation.cpp
  PLUGIN_XML
    DepthImageBasedGeometryApproximation.xml
  LINK
    depthImageBasedGeometryApproximation
    )

