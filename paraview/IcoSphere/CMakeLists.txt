ttk_add_paraview_plugin(ttkIcoSphere
  SOURCES
    ${VTKWRAPPER_DIR}/ttkIcoSphere/ttkIcoSphere.cpp
  PLUGIN_XML
    IcoSphere.xml
  LINK
    icoSphere
    )

