ttk_add_paraview_plugin(ttkTextureMapFromField
  SOURCES
    ${VTKWRAPPER_DIR}/ttkTextureMapFromField/ttkTextureMapFromField.cpp
  PLUGIN_XML
    TextureMapFromField.xml
  LINK
    common
    )

