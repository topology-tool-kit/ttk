# Allows to disable each filter
option(TTK_BUILD_FTRGRAPH_FILTER "Build the FTRGraph filter" ON)
mark_as_advanced(TTK_BUILD_FTRGRAPH_FILTER)

if(${TTK_BUILD_FTRGRAPH_FILTER})
  ttk_register_pv_filter(pvFTRGraph ttkFTRGraph)
endif()
