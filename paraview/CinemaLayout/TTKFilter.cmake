# Allows to disable each filter
option(TTK_BUILD_CINEMALAYOUT_FILTER "Build the CinemaLayout filter" ON)
mark_as_advanced(TTK_BUILD_CINEMALAYOUT_FILTER)

if(${TTK_BUILD_CINEMALAYOUT_FILTER})
  ttk_register_pv_filter(pvCinemaLayout ttkCinemaLayout)
endif()
