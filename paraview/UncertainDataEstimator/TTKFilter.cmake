# Allows to disable each filter
option(TTK_BUILD_UNCERTAINDATAESTIMATOR_FILTER "Build the UncertainDataEstimator filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_UNCERTAINDATAESTIMATOR_FILTER)

if(${TTK_BUILD_UNCERTAINDATAESTIMATOR_FILTER})
  ttk_register_pv_filter(ttkUncertainDataEstimator UncertainDataEstimator.xml)
endif()
