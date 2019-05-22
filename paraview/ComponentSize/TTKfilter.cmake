# Allows to disable each filter
option(TTK_BUILD_COMPONENTSIZE_FILTER "Build the ComponentSize filter" ON)
mark_as_advanced(TTK_BUILD_COMPONENTSIZE_FILTER)

if(${TTK_BUILD_COMPONENTSIZE_FILTER})
  ttk_register_pv_filter(pvComponentSize ttkComponentSize)
endif()
