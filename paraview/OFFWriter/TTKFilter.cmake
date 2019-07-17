# Allows to disable each filter
option(TTK_BUILD_OFFWRITER_FILTER "Build the OFFWriter filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_OFFWRITER_FILTER)

if(${TTK_BUILD_OFFWRITER_FILTER})
  ttk_register_pv_filter(pvOFFWriter ttkOFFWriter)
endif()
