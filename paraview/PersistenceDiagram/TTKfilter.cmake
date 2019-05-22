# Allows to disable each filter
option(TTK_BUILD_PERSISTENCEDIAGRAM_FILTER "Build the PersistenceDiagram filter" ON)
mark_as_advanced(TTK_BUILD_PERSISTENCEDIAGRAM_FILTER)

if(${TTK_BUILD_PERSISTENCEDIAGRAM_FILTER})
  ttk_register_pv_filter(pvPersistenceDiagram ttkPersistenceDiagram)
endif()
