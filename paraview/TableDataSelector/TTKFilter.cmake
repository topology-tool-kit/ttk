option(TTK_BUILD_TABLEDATASELECTOR_FILTER "Build the TableDataSelector filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_TABLEDATASELECTOR_FILTER)

if(${TTK_BUILD_TABLEDATASELECTOR_FILTER})
  ttk_register_pv_filter(ttkTableDataSelector TableDataSelector.xml)
endif()
