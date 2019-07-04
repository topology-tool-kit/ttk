ttk_add_paraview_plugin(ttkBarycentricSubdivision
  SOURCES
    ${VTKWRAPPER_DIR}/ttkBarycentricSubdivision/ttkBarycentricSubdivision.cpp
  PLUGIN_XML
    BarycentricSubdivision.xml
  LINK
    barycentricSubdivision
    )

