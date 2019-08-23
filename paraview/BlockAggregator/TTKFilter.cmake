option(TTK_BUILD_BLOCKAGGREGATOR_FILTER "Build the BlockAggregator filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_BLOCKAGGREGATOR_FILTER)

if(${TTK_BUILD_BLOCKAGGREGATOR_FILTER})
  ttk_register_pv_filter(ttkBlockAggregator BlockAggregator.xml)
endif()
