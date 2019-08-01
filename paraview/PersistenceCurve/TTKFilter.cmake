# Allows to disable each filter
option(TTK_BUILD_PERSISTENCECURVE_FILTER "Build the PersistenceCurve filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_PERSISTENCECURVE_FILTER)

if(${TTK_BUILD_PERSISTENCECURVE_FILTER})
  ttk_register_pv_filter(ttkPersistenceCurve PersistenceCurve.xml)
endif()
