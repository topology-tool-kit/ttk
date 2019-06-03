# Allows to disable each filter
option(TTK_BUILD_RANGEPOLYGON_FILTER "Build the RangePolygon filter" ON)
mark_as_advanced(TTK_BUILD_RANGEPOLYGON_FILTER)

if(${TTK_BUILD_RANGEPOLYGON_FILTER})
  ttk_register_pv_filter(pvRangePolygon ttkRangePolygon)
endif()
