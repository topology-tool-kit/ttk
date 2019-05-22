# Allows to disable each filter
option(TTK_BUILD_SPHEREFROMPOINT_FILTER "Build the SphereFromPoint filter" ON)
mark_as_advanced(TTK_BUILD_SPHEREFROMPOINT_FILTER)

if(${TTK_BUILD_SPHEREFROMPOINT_FILTER})
  ttk_register_pv_filter(pvSphereFromPoint ttkSphereFromPoint)
endif()
