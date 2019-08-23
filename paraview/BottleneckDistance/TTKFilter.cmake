option(TTK_BUILD_BOTTLENECKDISTANCE_FILTER "Build the BottleneckDistance filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_BOTTLENECKDISTANCE_FILTER)

if(${TTK_BUILD_BOTTLENECKDISTANCE_FILTER})
  ttk_register_pv_filter(ttkBottleneckDistance BottleneckDistance.xml)
endif()
