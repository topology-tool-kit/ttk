option(TTK_BUILD_BLANK_FILTER "Build the Blank filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_BLANK_FILTER)

if(${TTK_BUILD_BLANK_FILTER})
  # rmq, could use a ttk.plugin for these parameters
  ttk_register_pv_filter(ttkBlank Blank.xml)
endif()
