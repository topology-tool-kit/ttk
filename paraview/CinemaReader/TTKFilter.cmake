option(TTK_BUILD_CINEMAREADER_FILTER "Build the CinemaReader filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_CINEMAREADER_FILTER)

if(${TTK_BUILD_CINEMAREADER_FILTER})
  ttk_register_pv_filter(ttkCinemaReader CinemaReader.xml)
endif()
