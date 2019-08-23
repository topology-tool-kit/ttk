option(TTK_BUILD_CINEMAPRODUCTREADER_FILTER "Build the CinemaProductReader filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_CINEMAPRODUCTREADER_FILTER)

if(${TTK_BUILD_CINEMAPRODUCTREADER_FILTER})
  ttk_register_pv_filter(ttkCinemaProductReader CinemaProductReader.xml)
endif()
