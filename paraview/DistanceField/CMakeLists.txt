ttk_add_paraview_plugin(ttkDistanceField
  SOURCES
    ${VTKWRAPPER_DIR}/ttkDistanceField/ttkDistanceField.cpp
  PLUGIN_XML
    DistanceField.xml
  LINK
    distanceField
    )

