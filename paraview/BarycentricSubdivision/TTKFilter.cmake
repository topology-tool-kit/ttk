option(TTK_BUILD_BARYCENTRICSUBDIVISION_FILTER "Build the BarycentricSubdivision filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_BARYCENTRICSUBDIVISION_FILTER)

if(${TTK_BUILD_BARYCENTRICSUBDIVISION_FILTER})
  ttk_register_pv_filter(ttkBarycentricSubdivision BarycentricSubdivision.xml)
endif()

