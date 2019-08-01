# Allows to disable each filter
option(TTK_BUILD_PERSISTENCEDIAGRAM_FILTER "Build the PersistenceDiagram filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_PERSISTENCEDIAGRAM_FILTER)

if(${TTK_BUILD_PERSISTENCEDIAGRAM_FILTER})
  ttk_register_pv_filter(ttkPersistenceDiagram PersistenceDiagram.xml)
endif()
