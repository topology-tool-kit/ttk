ttk_add_paraview_plugin(ttkTrackingFromOverlap
	SOURCES ${VTKWRAPPER_DIR}/ttkTrackingFromOverlap/ttkTrackingFromOverlap.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML TrackingFromOverlap.xml
	LINK trackingFromOverlap)