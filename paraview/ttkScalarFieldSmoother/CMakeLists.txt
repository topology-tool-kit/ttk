ttk_add_paraview_plugin(ttkScalarFieldSmoother
	SOURCES ${VTKWRAPPER_DIR}/ttkScalarFieldSmoother/ttkScalarFieldSmoother.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML ScalarFieldSmoother.xml
	LINK scalarFieldSmoother ttkTriangulation)

