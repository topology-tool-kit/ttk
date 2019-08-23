option(TTK_BUILD_JACOBISET_FILTER "Build the JacobiSet filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_JACOBISET_FILTER)

if(${TTK_BUILD_JACOBISET_FILTER})
  ttk_register_pv_filter(ttkJacobiSet JacobiSet.xml)
endif()
