# Allows to disable each filter
option(TTK_BUILD_INTEGRALLINES_FILTER "Build the IntegralLines filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_INTEGRALLINES_FILTER)

if(${TTK_BUILD_INTEGRALLINES_FILTER})
  ttk_register_pv_filter(ttkIntegralLines IntegralLines.xml)
endif()
