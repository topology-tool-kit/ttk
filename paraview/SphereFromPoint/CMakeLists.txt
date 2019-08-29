ttk_add_paraview_plugin(ttkSphereFromPoint
  SOURCES
    ${VTKWRAPPER_DIR}/ttkSphereFromPoint/ttkSphereFromPoint.cpp
  PLUGIN_XML
    SphereFromPoint.xml
  LINK
    common
    )

