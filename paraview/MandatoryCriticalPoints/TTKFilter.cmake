# Allows to disable each filter
option(TTK_BUILD_MANDATORYCRITICALPOINTS_FILTER "Build the MandatoryCriticalPoints filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_MANDATORYCRITICALPOINTS_FILTER)

if(${TTK_BUILD_MANDATORYCRITICALPOINTS_FILTER})
  ttk_register_pv_filter(pvMandatoryCriticalPoints ttkMandatoryCriticalPoints)
endif()
