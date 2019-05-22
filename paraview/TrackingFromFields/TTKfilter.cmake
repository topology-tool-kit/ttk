# Allows to disable each filter
option(TTK_BUILD_TRACKINGFROMFIELDS_FILTER "Build the TrackingFromFields filter" ON)
mark_as_advanced(TTK_BUILD_TRACKINGFROMFIELDS_FILTER)

if(${TTK_BUILD_TRACKINGFROMFIELDS_FILTER})
  ttk_register_pv_filter(pvTrackingFromFields ttkTrackingFromFields)
endif()
