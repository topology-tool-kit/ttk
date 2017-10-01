ttk_add_paraview_plugin(ttkPersistenceDiagram
	SOURCES ${VTKWRAPPER_DIR}/ttkPersistenceDiagram/ttkPersistenceDiagram.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML PersistenceDiagram.xml
	LINK persistenceDiagram ttkTriangulation)

