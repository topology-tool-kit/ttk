option(TTK_BUILD_FTMTREE_FILTER "Build the FTMTree filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_FTMTREE_FILTER)

if(${TTK_BUILD_FTMTREE_FILTER})
  ttk_register_pv_filter(ttkFTMTree FTMTree.xml)
endif()
