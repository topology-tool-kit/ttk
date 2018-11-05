ttk_add_paraview_plugin(ttkCinemaQuery
	SOURCES ${VTKWRAPPER_DIR}/ttkCinemaQuery/ttkCinemaQuery.cpp
	PLUGIN_VERSION 1.0
	PLUGIN_XML CinemaQuery.xml
	LINK cinemaQuery)

