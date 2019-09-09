ttk_add_paraview_plugin(ttkImportEmbeddingFromTable
  SOURCES
    ${VTKWRAPPER_DIR}/ttkImportEmbeddingFromTable/ttkImportEmbeddingFromTable.cpp
  PLUGIN_XML
    ImportEmbeddingFromTable.xml
  LINK
    common
    )

