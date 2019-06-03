# Allows to disable each filter
option(TTK_BUILD_TOPOLOGICALCOMPRESSION_FILTER "Build the TopologicalCompression filter" ON)
mark_as_advanced(TTK_BUILD_TOPOLOGICALCOMPRESSION_FILTER)

if(${TTK_BUILD_TOPOLOGICALCOMPRESSION_FILTER})
  ttk_register_pv_filter(pvTopologicalCompression ttkTopologicalCompression)
endif()
