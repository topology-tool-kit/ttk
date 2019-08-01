# Allows to disable each filter
option(TTK_BUILD_DATASETTOTABLE_FILTER "Build the DataSetToTable filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_DATASETTOTABLE_FILTER)

if(${TTK_BUILD_DATASETTOTABLE_FILTER})
  ttk_register_pv_filter(ttkDataSetToTable DataSetToTable.xml)
endif()
