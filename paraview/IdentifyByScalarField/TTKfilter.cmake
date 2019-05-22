# Allows to disable each filter
option(TTK_BUILD_IDENTIFYBYSCALARFIELD_FILTER "Build the IdentifyByScalarField filter" ON)
mark_as_advanced(TTK_BUILD_IDENTIFYBYSCALARFIELD_FILTER)

if(${TTK_BUILD_IDENTIFYBYSCALARFIELD_FILTER})
  ttk_register_pv_filter(pvIdentifyByScalarField ttkIdentifyByScalarField)
endif()
