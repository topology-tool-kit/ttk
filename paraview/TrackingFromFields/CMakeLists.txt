include_directories(${VTKWRAPPER_DIR}/ttkBottleneckDistance/)
include_directories(${VTKWRAPPER_DIR}/ttkTrackingFromPersistenceDiagrams/)

ttk_add_paraview_plugin(ttkTrackingFromFields
  SOURCES
    ${VTKWRAPPER_DIR}/ttkTrackingFromFields/ttkTrackingFromFields.cpp
  PLUGIN_XML
    TrackingFromFields.xml
  LINK
    trackingFromFields
    trackingFromPersistenceDiagrams
    )
