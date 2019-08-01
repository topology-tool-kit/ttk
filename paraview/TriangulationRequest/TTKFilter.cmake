# Allows to disable each filter
option(TTK_BUILD_TRIANGULATIONREQUEST_FILTER "Build the TriangulationRequest filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_TRIANGULATIONREQUEST_FILTER)

if(${TTK_BUILD_TRIANGULATIONREQUEST_FILTER})
  ttk_register_pv_filter(ttkTriangulationRequest TriangulationRequest.xml)
endif()
