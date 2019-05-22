# Allows to disable each filter
option(TTK_BUILD_EIGENFIELD_FILTER "Build the EigenField filter" ON)
mark_as_advanced(TTK_BUILD_EIGENFIELD_FILTER)

if(${TTK_BUILD_EIGENFIELD_FILTER})
  ttk_register_pv_filter(pvEigenField ttkEigenField)
endif()
