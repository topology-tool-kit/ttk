# Allows to disable each filter
option(TTK_BUILD_QUADRANGULATION_SUBDIVISION_FILTER "Build the QuadrangulationSubdivision filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_QUADRANGULATION_SUBDIVISION_FILTER)

if(${TTK_BUILD_QUADRANGULATION_SUBDIVISION_FILTER})
  ttk_register_pv_filter(ttkQuadrangulationSubdivision QuadrangulationSubdivision.xml)
endif()
