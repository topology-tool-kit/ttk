# Allows to disable each filter
option(TTK_BUILD_TOPOLOGICALSIMPLIFICATION_FILTER "Build the TopologicalSimplification filter" ON)
mark_as_advanced(TTK_BUILD_TOPOLOGICALSIMPLIFICATION_FILTER)

if(${TTK_BUILD_TOPOLOGICALSIMPLIFICATION_FILTER})
  ttk_register_pv_filter(pvTopologicalSimplification ttkTopologicalSimplification)
endif()
