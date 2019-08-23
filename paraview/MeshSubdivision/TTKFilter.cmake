option(TTK_BUILD_MESHSUBDIVISION_FILTER "Build the MeshSubdivision filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_MESHSUBDIVISION_FILTER)

if(${TTK_BUILD_MESHSUBDIVISION_FILTER})
  ttk_register_pv_filter(ttkMeshSubdivision MeshSubdivision.xml)
endif()
