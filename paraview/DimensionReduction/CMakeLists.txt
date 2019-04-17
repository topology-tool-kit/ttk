ttk_add_paraview_plugin(ttkDimensionReduction
  SOURCES
    ${VTKWRAPPER_DIR}/ttkDimensionReduction/ttkDimensionReduction.cpp
  PLUGIN_XML
    DimensionReduction.xml
  LINK
    dimensionReduction
    )

