# Allows to disable each filter
option(TTK_BUILD_FOREACHROW_FILTER "Build the ForEachRow filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_FOREACHROW_FILTER)

if(${TTK_BUILD_FOREACHROW_FILTER})
  ttk_register_pv_filter(ttkForEachRow ForEachRow.xml)
endif()
