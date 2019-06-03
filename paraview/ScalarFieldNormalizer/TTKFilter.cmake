# Allows to disable each filter
option(TTK_BUILD_SCALARFIELDNORMALIZER_FILTER "Build the ScalarFieldNormalizer filter" ON)
mark_as_advanced(TTK_BUILD_SCALARFIELDNORMALIZER_FILTER)

if(${TTK_BUILD_SCALARFIELDNORMALIZER_FILTER})
  ttk_register_pv_filter(pvScalarFieldNormalizer ttkScalarFieldNormalizer)
endif()
