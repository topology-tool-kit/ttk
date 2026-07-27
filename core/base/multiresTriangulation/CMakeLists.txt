ttk_add_base_library(multiresTriangulation
  SOURCES
    MultiresTriangulation.cpp
    MultiresTriangulation_CustomDecimation.cpp
  HEADERS
    MultiresTriangulation.h
  DEPENDS
    triangulation
    )

