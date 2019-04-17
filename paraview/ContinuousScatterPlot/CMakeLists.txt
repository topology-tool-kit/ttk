ttk_add_paraview_plugin(ttkContinuousScatterPlot
  SOURCES
    ${VTKWRAPPER_DIR}/ttkContinuousScatterPlot/ttkContinuousScatterPlot.cpp
  PLUGIN_XML
    ContinuousScatterPlot.xml
  LINK
    continuousScatterPlot
    )
