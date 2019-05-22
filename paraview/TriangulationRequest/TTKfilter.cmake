# Allows to disable each filter
option(TTK_BUILD_TRIANGULATIONREQUEST_FILTER "Build the TriangulationRequest filter" ON)
mark_as_advanced(TTK_BUILD_TRIANGULATIONREQUEST_FILTER)

if(${TTK_BUILD_TRIANGULATIONREQUEST_FILTER})
  ttk_register_pv_filter(pvTriangulationRequest ttkTriangulationRequest)
endif()
