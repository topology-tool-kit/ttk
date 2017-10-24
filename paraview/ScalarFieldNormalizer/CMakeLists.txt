ttk_add_paraview_plugin(ttkScalarFieldNormalizer
	SOURCES ${VTKWRAPPER_DIR}/ttkScalarFieldNormalizer/ttkScalarFieldNormalizer.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML ScalarFieldNormalizer.xml
	LINK common)

