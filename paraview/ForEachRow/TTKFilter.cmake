# Allows to disable each filter
option(TTK_BUILD_FOREACHROW_FILTER "Build the ForEachRow filter" ON)
mark_as_advanced(TTK_BUILD_FOREACHROW_FILTER)

if(${TTK_BUILD_FOREACHROW_FILTER})
  ttk_register_pv_filter(pvForEachRow ttkForEachRow)
endif()
