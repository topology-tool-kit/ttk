option(TTK_BUILD_OFFREADER_FILTER "Build the OFFReader filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_OFFREADER_FILTER)

if(${TTK_BUILD_OFFREADER_FILTER})
  ttk_register_pv_filter(ttkOFFReader OFFReader.xml)
endif()
