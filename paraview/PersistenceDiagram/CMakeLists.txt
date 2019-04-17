ttk_add_paraview_plugin(ttkPersistenceDiagram
  SOURCES
    ${VTKWRAPPER_DIR}/ttkPersistenceDiagram/ttkPersistenceDiagram.cpp
  PLUGIN_XML
    PersistenceDiagram.xml
  LINK
    persistenceDiagram
    )

