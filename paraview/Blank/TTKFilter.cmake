# Allows to disable each filter
option(TTK_BUILD_BLANK_FILTER "Build the Blank filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_BLANK_FILTER)

if(${TTK_BUILD_BLANK_FILTER})
  ttk_register_pv_filter(pvBlank ttkBlank)
endif()

