# Allows to disable each filter
option(TTK_BUILD_POINTDATASELECTOR_FILTER "Build the PointDataSelector filter" ON)
mark_as_advanced(TTK_BUILD_POINTDATASELECTOR_FILTER)

if(${TTK_BUILD_POINTDATASELECTOR_FILTER})
  ttk_register_pv_filter(pvPointDataSelector ttkPointDataSelector)
endif()
