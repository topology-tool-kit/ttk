ttk_add_paraview_plugin(ttkGeometrySmoother
  SOURCES
    ${VTKWRAPPER_DIR}/ttkGeometrySmoother/ttkGeometrySmoother.cpp
  PLUGIN_XML
    GeometrySmoother.xml
  LINK
    scalarFieldSmoother
    )

