# Allows to disable each filter
option(TTK_BUILD_ICOSPHERE_FILTER "Build the IcoSphere filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_ICOSPHERE_FILTER)

if(${TTK_BUILD_ICOSPHERE_FILTER})
  ttk_register_pv_filter(ttkIcoSphere IcoSphere.xml)
endif()
