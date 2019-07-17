# Allows to disable each filter
option(TTK_BUILD_IDENTIFIERS_FILTER "Build the Identifiers filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_IDENTIFIERS_FILTER)

if(${TTK_BUILD_IDENTIFIERS_FILTER})
  ttk_register_pv_filter(pvIdentifiers ttkIdentifiers)
endif()
