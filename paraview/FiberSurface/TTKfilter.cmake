# Allows to disable each filter
option(TTK_BUILD_FIBERSURFACE_FILTER "Build the FiberSurface filter" ON)
mark_as_advanced(TTK_BUILD_FIBERSURFACE_FILTER)

if(${TTK_BUILD_FIBERSURFACE_FILTER})
  ttk_register_pv_filter(pvFiberSurface ttkFiberSurface)
endif()
