# Allows to disable each filter
option(TTK_BUILD_LDISTANCE_FILTER "Build the LDistance filter" ON)
mark_as_advanced(TTK_BUILD_LDISTANCE_FILTER)

if(${TTK_BUILD_LDISTANCE_FILTER})
  ttk_register_pv_filter(pvLDistance ttkLDistance)
endif()
