option(TTK_BUILD_DATASETINTERPOLATOR_FILTER "Build the DataSetInterpolator filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_DATASETINTERPOLATOR_FILTER)

if(${TTK_BUILD_DATASETINTERPOLATOR_FILTER})
  ttk_register_pv_filter(ttkDataSetInterpolator DataSetInterpolator.xml)
endif()
