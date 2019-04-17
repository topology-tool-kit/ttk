ttk_add_paraview_plugin(ttkMeshSubdivision
  SOURCES
    ${VTKWRAPPER_DIR}/ttkMeshSubdivision/ttkMeshSubdivision.cpp
  PLUGIN_XML
    MeshSubdivision.xml
  LINK
    common
    )

