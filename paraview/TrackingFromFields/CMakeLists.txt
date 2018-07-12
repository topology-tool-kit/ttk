ttk_add_paraview_plugin(ttkTrackingFromFields
  SOURCES ${VTKWRAPPER_DIR}/ttkTrackingFromFields/ttkTrackingFromFields.cpp
  PLUGIN_VERSION 1.0
  PLUGIN_XML TrackingFromFields.xml
  LINK ttkPersistenceDiagram ttkBottleneckDistance ttkTrackingFromPersistenceDiagrams ttkFTMTree)
