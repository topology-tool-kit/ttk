ttk_add_paraview_plugin(ttkUncertainDataEstimator
  SOURCES
    ${VTKWRAPPER_DIR}/ttkUncertainDataEstimator/ttkUncertainDataEstimator.cpp
  PLUGIN_XML
    UncertainDataEstimator.xml
  LINK
    uncertainDataEstimator
    )

