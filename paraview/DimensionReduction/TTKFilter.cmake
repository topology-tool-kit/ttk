option(TTK_BUILD_DIMENSIONREDUCTION_FILTER "Build the DimensionReduction filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_DIMENSIONREDUCTION_FILTER)

if(${TTK_BUILD_DIMENSIONREDUCTION_FILTER})
  ttk_register_pv_filter(ttkDimensionReduction DimensionReduction.xml)
endif()
