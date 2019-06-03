# Allows to disable each filter
option(TTK_BUILD_ENDFOR_FILTER "Build the EndFor filter" ON)
mark_as_advanced(TTK_BUILD_ENDFOR_FILTER)

if(${TTK_BUILD_ENDFOR_FILTER})
  ttk_register_pv_filter(pvEndFor ttkEndFor)
endif()
