# Allows to disable each filter
option(TTK_BUILD_CELLDATASELECTOR_FILTER "Build the CellDataSelector filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_CELLDATASELECTOR_FILTER)

if(${TTK_BUILD_CELLDATASELECTOR_FILTER})
  ttk_register_pv_filter(pvCellDataSelector ttkCellDataSelector)
endif()
