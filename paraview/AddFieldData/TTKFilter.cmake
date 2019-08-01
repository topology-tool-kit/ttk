# Allows to disable each filter
option(TTK_BUILD_ADDFIELDDATA_FILTER "Build the AddFieldData filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_ADDFIELDDATA_FILTER)

if(${TTK_BUILD_ADDFIELDDATA_FILTER})
  ttk_register_pv_filter(ttkAddFieldData AddFieldData.xml)
endif()
