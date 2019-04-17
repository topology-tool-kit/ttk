ttk_add_paraview_plugin(ttkJacobiSet
  SOURCES
    ${VTKWRAPPER_DIR}/ttkJacobiSet/ttkJacobiSet.cpp
  PLUGIN_XML
    JacobiSet.xml
  LINK
    jacobiSet
    )

