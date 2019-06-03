# Allows to disable each filter
option(TTK_BUILD_DISCRETEGRADIENT_FILTER "Build the DiscreteGradient filter" ON)
mark_as_advanced(TTK_BUILD_DISCRETEGRADIENT_FILTER)

if(${TTK_BUILD_DISCRETEGRADIENT_FILTER})
  ttk_register_pv_filter(pvDiscreteGradient ttkDiscreteGradient)
endif()
