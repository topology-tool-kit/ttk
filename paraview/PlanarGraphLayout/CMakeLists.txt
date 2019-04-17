ttk_add_paraview_plugin(ttkPlanarGraphLayout
  SOURCES
    ${VTKWRAPPER_DIR}/ttkPlanarGraphLayout/ttkPlanarGraphLayout.cpp
  PLUGIN_XML
    PlanarGraphLayout.xml
  LINK
    planarGraphLayout
    )

