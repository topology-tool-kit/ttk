ttk_add_paraview_plugin(ttkDiscreteGradient
	SOURCES ${VTKWRAPPER_DIR}/ttkDiscreteGradient/ttkDiscreteGradient.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML DiscreteGradient.xml
	LINK discreteGradient)

