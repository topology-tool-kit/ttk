ttk_add_paraview_plugin(ttkCinemaQuery
  SOURCES
    ${VTKWRAPPER_DIR}/ttkCinemaQuery/ttkCinemaQuery.cpp
  PLUGIN_XML
    CinemaQuery.xml
  LINK
    cinemaQuery
    )

