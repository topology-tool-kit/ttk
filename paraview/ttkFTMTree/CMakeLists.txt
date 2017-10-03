ttk_add_paraview_plugin(ttkFTMTree
	SOURCES ${VTKWRAPPER_DIR}/ttkFTMTree/ttkFTMTree.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML FTMTree.xml
	LINK ftmTree ttkTriangulation)

