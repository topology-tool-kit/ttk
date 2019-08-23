option(TTK_BUILD_CONTOURFORESTS_FILTER "Build the ContourForests filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_CONTOURFORESTS_FILTER)

if(${TTK_BUILD_CONTOURFORESTS_FILTER})
  ttk_register_pv_filter(ttkContourForests ContourForests.xml)
endif()
