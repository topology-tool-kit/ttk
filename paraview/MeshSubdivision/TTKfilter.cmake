# Allows to disable each filter
option(TTK_BUILD_MESHSUBDIVISION_FILTER "Build the MeshSubdivision filter" ON)
mark_as_advanced(TTK_BUILD_MESHSUBDIVISION_FILTER)

if(${TTK_BUILD_MESHSUBDIVISION_FILTER})
  ttk_register_pv_filter(pvMeshSubdivision ttkMeshSubdivision)
endif()
