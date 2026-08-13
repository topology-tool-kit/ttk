ttk_add_base_library(cinemaImaging
  SOURCES
    CinemaImagingEmbree.cpp
  HEADERS
    CinemaImaging.h
    CinemaImagingEmbree.h
    CinemaImagingNative.h
  DEPENDS
    boundingVolumeHierarchy
    triangulation
)

if (TTK_ENABLE_EMBREE AND EMBREE_FOUND)
  target_compile_definitions(cinemaImaging PUBLIC TTK_ENABLE_EMBREE)
  if(EMBREE_VERSION_MAJOR EQUAL 4)
    target_compile_definitions(cinemaImaging PUBLIC TTK_EMBREE4)
  endif()
  target_include_directories(cinemaImaging PUBLIC ${EMBREE_INCLUDE_DIRS})
  target_link_libraries(cinemaImaging PUBLIC ${EMBREE_LIBRARY})
endif()
