# Allows to disable each filter
option(TTK_BUILD_GAUSSIANPOINTCLOUD_FILTER "Build the GaussianPointCloud filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_GAUSSIANPOINTCLOUD_FILTER)

if(${TTK_BUILD_GAUSSIANPOINTCLOUD_FILTER})
  ttk_register_pv_filter(ttkGaussianPointCloud GaussianPointCloud.xml)
endif()
