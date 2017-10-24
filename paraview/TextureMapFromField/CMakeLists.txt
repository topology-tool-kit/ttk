ttk_add_paraview_plugin(ttkTextureMapFromField
	SOURCES ${VTKWRAPPER_DIR}/ttkTextureMapFromField/ttkTextureMapFromField.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML TextureMapFromField.xml
	LINK common)

