# Allows to disable each filter
option(TTK_BUILD_CINEMAWRITER_FILTER "Build the CinemaWriter filter" ON)
mark_as_advanced(TTK_BUILD_CINEMAWRITER_FILTER)

if(${TTK_BUILD_CINEMAWRITER_FILTER})
  ttk_register_pv_filter(pvCinemaWriter ttkCinemaWriter)
endif()
