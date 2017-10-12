ttk_add_paraview_plugin(ttkGeometrySmoother
	SOURCES ${VTKWRAPPER_DIR}/ttkGeometrySmoother/ttkGeometrySmoother.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML GeometrySmoother.xml
	LINK scalarFieldSmoother)

