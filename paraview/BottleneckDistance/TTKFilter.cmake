# Allows to disable each filter
option(TTK_BUILD_BOTTLENECKDISTANCE_FILTER "Build the BottleneckDistance filter" ON)
mark_as_advanced(TTK_BUILD_BOTTLENECKDISTANCE_FILTER)

if(${TTK_BUILD_BOTTLENECKDISTANCE_FILTER})
  ttk_register_pv_filter(pvBottleneckDistance ttkBottleneckDistance)
endif()
