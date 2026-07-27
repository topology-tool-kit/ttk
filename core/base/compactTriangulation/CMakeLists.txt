
ttk_add_base_library(compactTriangulation
  SOURCES
    CompactTriangulation.cpp
  HEADERS
    CompactTriangulation.h
  DEPENDS
    abstractTriangulation
    skeleton
    Boost::boost
    )
