# Allows to disable each filter
option(TTK_BUILD_POINTDATACONVERTER_FILTER "Build the PointDataConverter filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_POINTDATACONVERTER_FILTER)

if(${TTK_BUILD_POINTDATACONVERTER_FILTER})
  ttk_register_pv_filter(ttkPointDataConverter PointDataConverter.xml)
endif()
