# Allows to disable each filter
option(TTK_BUILD_REEBSPACE_FILTER "Build the ReebSpace filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_REEBSPACE_FILTER)

if(${TTK_BUILD_REEBSPACE_FILTER})
  ttk_register_pv_filter(ttkReebSpace ReebSpace.xml)
endif()
