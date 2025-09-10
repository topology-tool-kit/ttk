ttk_add_base_library(delaunayRipsPersistenceDiagram
  SOURCES
    DelaunayRipsPersistenceDiagram.cpp
  HEADERS
    DelaunayRipsPersistenceDiagram.h
    geoPH3.h
    geoPHd.h
    geoPHUtils.h
  DEPENDS
    ripsPersistenceDiagram
)

set(TTK_DELAUNAY_MAXIMUM_DIMENSION "8" CACHE STRING "Maximum compiled dimension for arbitrary dimension Delaunay-related code")
mark_as_advanced(TTK_DELAUNAY_MAXIMUM_DIMENSION)
target_compile_definitions(delaunayRipsPersistenceDiagram PUBLIC TTK_DELAUNAY_MAXIMUM_DIMENSION=${TTK_DELAUNAY_MAXIMUM_DIMENSION})

target_link_libraries(delaunayRipsPersistenceDiagram PUBLIC Boost::boost)

if(TTK_ENABLE_EIGEN)
  target_compile_definitions(delaunayRipsPersistenceDiagram PUBLIC TTK_ENABLE_EIGEN)
  target_include_directories(delaunayRipsPersistenceDiagram PUBLIC ${EIGEN3_INCLUDE_DIR})
  target_link_libraries(delaunayRipsPersistenceDiagram PUBLIC ${EIGEN3_LIBRARIES})
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
  target_link_libraries (delaunayRipsPersistenceDiagram PUBLIC TBB::tbb)
  add_compile_definitions(ENABLE_TBB)
  add_compile_definitions(CGAL_LINKED_WITH_TBB)
endif()