option(TTK_BUILD_RANGEPOLYGON_FILTER "Build the RangePolygon filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_RANGEPOLYGON_FILTER)

if(${TTK_BUILD_RANGEPOLYGON_FILTER})
  ttk_register_pv_filter(ttkRangePolygon RangePolygon.xml)
endif()
