ttk_add_base_library(delaunayRipsPersistenceDiagram
  SOURCES
    DelaunayRipsPersistenceDiagram.cpp
  HEADERS
    DelaunayRipsPersistenceDiagram.h
    dset.h
    geoPH3.h
    geoPHd.h
    geoPHUtils.h
  DEPENDS
    ripsPersistenceDiagram
)

set(TTK_DELAUNAY_MAX_COMPILED_DIMENSION "8" CACHE STRING "Maximum compiled dimension for arbitrary dimension Delaunay-related code")
mark_as_advanced(TTK_DELAUNAY_MAX_COMPILED_DIMENSION)
target_compile_definitions(delaunayRipsPersistenceDiagram PUBLIC TTK_DELAUNAY_MAX_COMPILED_DIMENSION=${TTK_DELAUNAY_MAX_COMPILED_DIMENSION})

target_link_libraries(delaunayRipsPersistenceDiagram PUBLIC Boost::boost)

if(TTK_ENABLE_EIGEN)
  target_compile_definitions(delaunayRipsPersistenceDiagram PUBLIC TTK_ENABLE_EIGEN)
  target_link_libraries(delaunayRipsPersistenceDiagram PUBLIC Eigen3::Eigen)
endif()

if(TTK_ENABLE_CGAL)
  target_compile_definitions(delaunayRipsPersistenceDiagram PUBLIC TTK_ENABLE_CGAL)
  if(CGAL_VERSION VERSION_LESS "6.0")
    target_link_libraries(delaunayRipsPersistenceDiagram PUBLIC ${CGAL_LIBRARIES})
  else()
    target_link_libraries(delaunayRipsPersistenceDiagram PUBLIC CGAL::CGAL)
  endif()
endif()

if(TTK_ENABLE_TBB)
  target_compile_definitions(delaunayRipsPersistenceDiagram PUBLIC TTK_ENABLE_TBB)
  target_compile_definitions(delaunayRipsPersistenceDiagram PUBLIC CGAL_LINKED_WITH_TBB)
  target_link_libraries(delaunayRipsPersistenceDiagram PUBLIC TBB::tbb TBB::tbbmalloc)
endif()