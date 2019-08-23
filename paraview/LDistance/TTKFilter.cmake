option(TTK_BUILD_LDISTANCE_FILTER "Build the LDistance filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_LDISTANCE_FILTER)

if(${TTK_BUILD_LDISTANCE_FILTER})
  ttk_register_pv_filter(ttkLDistance LDistance.xml)
endif()
