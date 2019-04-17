ttk_add_paraview_plugin(ttkFiberSurface
  SOURCES
    ${VTKWRAPPER_DIR}/ttkFiberSurface/ttkFiberSurface.cpp
  PLUGIN_XML
    FiberSurface.xml
  LINK
    fiberSurface
    )


