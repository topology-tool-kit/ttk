# Allows to disable each filter
option(TTK_BUILD_OBJWRITER_FILTER "Build the OBJWriter filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_OBJWRITER_FILTER)

if(${TTK_BUILD_OBJWRITER_FILTER})
  ttk_register_pv_filter(pvOBJWriter ttkOBJWriter)
endif()
