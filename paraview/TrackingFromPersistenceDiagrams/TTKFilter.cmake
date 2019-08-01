# Allows to disable each filter
option(TTK_BUILD_TRACKINGFROMPERSISTENCEDIAGRAMS_FILTER "Build the TrackingFromPersistenceDiagrams filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_TRACKINGFROMPERSISTENCEDIAGRAMS_FILTER)

if(${TTK_BUILD_TRACKINGFROMPERSISTENCEDIAGRAMS_FILTER})
  ttk_register_pv_filter(ttkTrackingFromPersistenceDiagrams TrackingFromPersistenceDiagrams.xml)
endif()
