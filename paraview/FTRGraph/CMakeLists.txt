ttk_add_paraview_plugin(ttkFTRGraph
	SOURCES ${VTKWRAPPER_DIR}/ttkFTRGraph/ttkFTRGraph.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML FTRGraph.xml
	LINK ftrGraph)

