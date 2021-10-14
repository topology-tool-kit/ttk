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
  target_include_directories(cinemaImaging PUBLIC ${EMBREE_INCLUDE_DIR})
  target_link_libraries(cinemaImaging PUBLIC ${EMBREE_LIBRARY})
endif()
