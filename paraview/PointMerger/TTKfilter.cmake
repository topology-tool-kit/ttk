# Allows to disable each filter
option(TTK_BUILD_POINTMERGER_FILTER "Build the PointMerger filter" ON)
mark_as_advanced(TTK_BUILD_POINTMERGER_FILTER)

if(${TTK_BUILD_POINTMERGER_FILTER})
  ttk_register_pv_filter(pvPointMerger ttkPointMerger)
endif()
