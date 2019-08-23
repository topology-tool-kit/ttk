option(TTK_BUILD_FIBERSURFACE_FILTER "Build the FiberSurface filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_FIBERSURFACE_FILTER)

if(${TTK_BUILD_FIBERSURFACE_FILTER})
  ttk_register_pv_filter(ttkFiberSurface FiberSurface.xml)
endif()
