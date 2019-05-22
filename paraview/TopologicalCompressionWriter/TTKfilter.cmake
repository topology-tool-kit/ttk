# Allows to disable each filter
option(TTK_BUILD_TOPOLOGICALCOMPRESSIONWRITER_FILTER "Build the TopologicalCompressionWriter filter" ON)
mark_as_advanced(TTK_BUILD_TOPOLOGICALCOMPRESSIONWRITER_FILTER)

if(${TTK_BUILD_TOPOLOGICALCOMPRESSIONWRITER_FILTER})
  ttk_register_pv_filter(pvTopologicalCompressionWriter ttkTopologicalCompressionWriter)
endif()
