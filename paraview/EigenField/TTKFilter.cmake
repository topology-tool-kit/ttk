# Allows to disable each filter
option(TTK_BUILD_EIGENFIELD_FILTER "Build the EigenField filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_EIGENFIELD_FILTER)

if(${TTK_BUILD_EIGENFIELD_FILTER})
  ttk_register_pv_filter(ttkEigenField EigenField.xml)
endif()
