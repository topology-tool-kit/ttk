ttk_add_paraview_plugin(ttkDiscreteGradient
  SOURCES
    ${VTKWRAPPER_DIR}/ttkDiscreteGradient/ttkDiscreteGradient.cpp
  PLUGIN_XML
    DiscreteGradient.xml
  LINK
    discreteGradient
    )
