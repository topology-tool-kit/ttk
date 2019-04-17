ttk_add_paraview_plugin(ttkTrackingFromOverlap
  SOURCES
    ${VTKWRAPPER_DIR}/ttkTrackingFromOverlap/ttkTrackingFromOverlap.cpp
  PLUGIN_XML
    TrackingFromOverlap.xml
  LINK
    trackingFromOverlap
    )
