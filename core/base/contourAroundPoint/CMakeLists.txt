# ttk_add_baseCode_package(unionFind)
ttk_add_base_library(contourAroundPoint
  SOURCES ContourAroundPoint.cpp
  HEADERS ContourAroundPoint.hpp
  DEPENDS triangulation)

# set_property(TARGET contourAroundPoint PROPERTY CXX_STANDARD 14)
