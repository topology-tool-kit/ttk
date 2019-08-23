option(TTK_BUILD_CINEMALAYOUT_FILTER "Build the CinemaLayout filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_CINEMALAYOUT_FILTER)

if(${TTK_BUILD_CINEMALAYOUT_FILTER})
  ttk_register_pv_filter(ttkCinemaLayout CinemaLayout.xml)
endif()
