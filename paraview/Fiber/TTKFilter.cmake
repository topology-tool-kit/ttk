option(TTK_BUILD_FIBER_FILTER "Build the Fiber filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_FIBER_FILTER)

if(${TTK_BUILD_FIBER_FILTER})
  ttk_register_pv_filter(ttkFiber Fiber.xml)
endif()
