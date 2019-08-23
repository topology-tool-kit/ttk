option(TTK_BUILD_ENDFOR_FILTER "Build the EndFor filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_ENDFOR_FILTER)

if(${TTK_BUILD_ENDFOR_FILTER})
  ttk_register_pv_filter(ttkEndFor EndFor.xml)
endif()
