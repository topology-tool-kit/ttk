option(TTK_BUILD_MORSESMALECOMPLEX_FILTER "Build the MorseSmaleComplex filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_MORSESMALECOMPLEX_FILTER)

if(${TTK_BUILD_MORSESMALECOMPLEX_FILTER})
  ttk_register_pv_filter(ttkMorseSmaleComplex MorseSmaleComplex.xml)
endif()
