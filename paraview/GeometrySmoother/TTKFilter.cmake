option(TTK_BUILD_GEOMETRYSMOOTHER_FILTER "Build the GeometrySmoother filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_GEOMETRYSMOOTHER_FILTER)

if(${TTK_BUILD_GEOMETRYSMOOTHER_FILTER})
  ttk_register_pv_filter(ttkGeometrySmoother GeometrySmoother.xml)
endif()
