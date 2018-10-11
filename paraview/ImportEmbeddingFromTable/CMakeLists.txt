ttk_add_paraview_plugin(ttkImportEmbeddingFromTable
	SOURCES ${VTKWRAPPER_DIR}/ttkImportEmbeddingFromTable/ttkImportEmbeddingFromTable.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML ImportEmbeddingFromTable.xml
	LINK common)

