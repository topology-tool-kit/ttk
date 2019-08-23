option(TTK_BUILD_POINTMERGER_FILTER "Build the PointMerger filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_POINTMERGER_FILTER)

if(${TTK_BUILD_POINTMERGER_FILTER})
  ttk_register_pv_filter(ttkPointMerger PointMerger.xml)
endif()
