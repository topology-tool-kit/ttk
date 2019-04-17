ttk_add_paraview_plugin(ttkMeshGraph
  SOURCES
    ${VTKWRAPPER_DIR}/ttkMeshGraph/ttkMeshGraph.cpp
  PLUGIN_XML
    MeshGraph.xml
  LINK
    meshGraph
    )

