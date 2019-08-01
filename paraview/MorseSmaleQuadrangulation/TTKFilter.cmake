# Allows to disable each filter
option(TTK_BUILD_MORSE_SMALE_QUADRANGULATION_FILTER "Build the MorseSmaleQuadrangulation filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_MORSE_SMALE_QUADRANGULATION_FILTER)

if(${TTK_BUILD_MORSE_SMALE_QUADRANGULATION_FILTER})
  ttk_register_pv_filter(ttkMorseSmaleQuadrangulation MorseSmaleQuadrangulation.xml)
endif()
