# Allows to disable each filter
option(TTK_BUILD_MANIFOLDCHECK_FILTER "Build the ManifoldCheck filter" ON)
mark_as_advanced(TTK_BUILD_MANIFOLDCHECK_FILTER)

if(${TTK_BUILD_MANIFOLDCHECK_FILTER})
  ttk_register_pv_filter(pvManifoldCheck ttkManifoldCheck)
endif()
