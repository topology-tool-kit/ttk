ttk_add_paraview_plugin(ttkMandatoryCriticalPoints
  SOURCES
    ${VTKWRAPPER_DIR}/ttkMandatoryCriticalPoints/ttkMandatoryCriticalPoints.cpp
  PLUGIN_XML
    MandatoryCriticalPoints.xml
  LINK
    mandatoryCriticalPoints
    )

