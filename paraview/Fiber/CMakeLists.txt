ttk_add_paraview_plugin(ttkFiber
  SOURCES
    ${VTKWRAPPER_DIR}/ttkFiber/ttkFiber.cpp
  PLUGIN_XML
    Fiber.xml
  LINK
    common
    )

