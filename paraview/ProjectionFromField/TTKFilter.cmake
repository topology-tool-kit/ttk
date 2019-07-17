# Allows to disable each filter
option(TTK_BUILD_PROJECTIONFROMFIELD_FILTER "Build the ProjectionFromField filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_PROJECTIONFROMFIELD_FILTER)

if(${TTK_BUILD_PROJECTIONFROMFIELD_FILTER})
  ttk_register_pv_filter(pvProjectionFromField ttkProjectionFromField)
endif()
