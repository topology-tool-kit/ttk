option(TTK_BUILD_TRIANGULATIONFILTER_FILTER "Build the TriangulationFilter filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_TRIANGULATIONFILTER_FILTER)

if(${TTK_BUILD_TRIANGULATIONFILTER_FILTER})
  ttk_register_pv_module(ttkTriangulationFilter)
endif()
