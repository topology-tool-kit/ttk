# Allows to disable each filter
option(TTK_BUILD_SCALARFIELDSMOOTHER_FILTER "Build the ScalarFieldSmoother filter" ON)
mark_as_advanced(TTK_BUILD_SCALARFIELDSMOOTHER_FILTER)

if(${TTK_BUILD_SCALARFIELDSMOOTHER_FILTER})
  ttk_register_pv_filter(pvScalarFieldSmoother ttkScalarFieldSmoother)
endif()
