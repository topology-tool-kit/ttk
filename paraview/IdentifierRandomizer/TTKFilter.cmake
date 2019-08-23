option(TTK_BUILD_IDENTIFIERRANDOMIZER_FILTER "Build the IdentifierRandomizer filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_IDENTIFIERRANDOMIZER_FILTER)

if(${TTK_BUILD_IDENTIFIERRANDOMIZER_FILTER})
  ttk_register_pv_filter(ttkIdentifierRandomizer IdentifierRandomizer.xml)
endif()
