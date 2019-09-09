ttk_add_paraview_plugin(ttkTopologicalSimplification
  SOURCES
    ${VTKWRAPPER_DIR}/ttkTopologicalSimplification/ttkTopologicalSimplification.cpp
  PLUGIN_XML
    TopologicalSimplification.xml
  LINK
    topologicalSimplification
    )

