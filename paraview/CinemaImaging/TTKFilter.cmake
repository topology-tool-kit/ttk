option(TTK_BUILD_CINEMAIMAGING_FILTER "Build the CinemaImaging filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_CINEMAIMAGING_FILTER)

if(${TTK_BUILD_CINEMAIMAGING_FILTER})
  ttk_register_pv_filter(ttkCinemaImaging CinemaImaging.xml)
endif()
