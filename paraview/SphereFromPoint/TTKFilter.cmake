# Allows to disable each filter
option(TTK_BUILD_SPHEREFROMPOINT_FILTER "Build the SphereFromPoint filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_SPHEREFROMPOINT_FILTER)

if(${TTK_BUILD_SPHEREFROMPOINT_FILTER})
  ttk_register_pv_filter(ttkSphereFromPoint SphereFromPoint.xml)
endif()
