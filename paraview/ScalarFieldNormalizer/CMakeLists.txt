ttk_add_paraview_plugin(ttkScalarFieldNormalizer
  SOURCES
    ${VTKWRAPPER_DIR}/ttkScalarFieldNormalizer/ttkScalarFieldNormalizer.cpp
  PLUGIN_XML
    ScalarFieldNormalizer.xml
  LINK
    common
    )

