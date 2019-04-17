ttk_add_paraview_plugin(ttkGaussianPointCloud
  SOURCES
    ${VTKWRAPPER_DIR}/ttkGaussianPointCloud/ttkGaussianPointCloud.cpp
  PLUGIN_XML
    GaussianPointCloud.xml
  LINK
    gaussianPointCloud
    )

