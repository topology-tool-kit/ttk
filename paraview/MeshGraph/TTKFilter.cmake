option(TTK_BUILD_MESHGRAPH_FILTER "Build the MeshGraph filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_MESHGRAPH_FILTER)

if(${TTK_BUILD_MESHGRAPH_FILTER})
  ttk_register_pv_filter(ttkMeshGraph MeshGraph.xml)
endif()
