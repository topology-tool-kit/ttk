include_directories(${VTKWRAPPER_DIR}/ttkBottleneckDistance/)

ttk_add_paraview_plugin(ttkTrackingFromPersistenceDiagrams
  SOURCES
    ${VTKWRAPPER_DIR}/ttkTrackingFromPersistenceDiagrams/ttkTrackingFromPersistenceDiagrams.cpp
  PLUGIN_XML
    TrackingFromPersistenceDiagrams.xml
  LINK
    trackingFromPersistenceDiagrams
    )

