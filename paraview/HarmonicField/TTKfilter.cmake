# Allows to disable each filter
option(TTK_BUILD_HARMONICFIELD_FILTER "Build the HarmonicField filter" ON)
mark_as_advanced(TTK_BUILD_HARMONICFIELD_FILTER)

if(${TTK_BUILD_HARMONICFIELD_FILTER})
  ttk_register_pv_filter(pvHarmonicField ttkHarmonicField)
endif()
