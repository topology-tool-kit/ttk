option(TTK_BUILD_CINEMAQUERY_FILTER "Build the CinemaQuery filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_CINEMAQUERY_FILTER)

if(${TTK_BUILD_CINEMAQUERY_FILTER})
  ttk_register_pv_filter(ttkCinemaQuery CinemaQuery.xml)
endif()
