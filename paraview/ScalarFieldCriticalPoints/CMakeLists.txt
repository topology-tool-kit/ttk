ttk_add_paraview_plugin(ttkScalarFieldCriticalPoints
  SOURCES
    ${VTKWRAPPER_DIR}/ttkScalarFieldCriticalPoints/ttkScalarFieldCriticalPoints.cpp
  PLUGIN_XML
    ScalarFieldCriticalPoints.xml
  LINK
    scalarFieldCriticalPoints
    )

