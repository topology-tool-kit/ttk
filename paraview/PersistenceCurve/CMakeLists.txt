ttk_add_paraview_plugin(ttkPersistenceCurve
  SOURCES
    ${VTKWRAPPER_DIR}/ttkPersistenceCurve/ttkPersistenceCurve.cpp
  PLUGIN_XML
    PersistenceCurve.xml
  LINK
    persistenceCurve
    )

