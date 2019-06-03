# Allows to disable each filter
option(TTK_BUILD_TABLEDATASELECTOR_FILTER "Build the TableDataSelector filter" ON)
mark_as_advanced(TTK_BUILD_TABLEDATASELECTOR_FILTER)

if(${TTK_BUILD_TABLEDATASELECTOR_FILTER})
  ttk_register_pv_filter(pvTableDataSelector ttkTableDataSelector)
endif()
