# Allows to disable each filter
option(TTK_BUILD_MANDATORYCRITICALPOINTS_FILTER "Build the MandatoryCriticalPoints filter" ON)
mark_as_advanced(TTK_BUILD_MANDATORYCRITICALPOINTS_FILTER)

if(${TTK_BUILD_MANDATORYCRITICALPOINTS_FILTER})
  ttk_register_pv_filter(pvMandatoryCriticalPoints ttkMandatoryCriticalPoints)
endif()
