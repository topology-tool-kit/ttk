ttk_add_base_library(triangulation
  SOURCES
    Triangulation.cpp
  HEADERS
    Triangulation.h
  DEPENDS
    abstractTriangulation
    compactTriangulation
    explicitTriangulation
    regularGridTriangulation    
    implicitTriangulation
    periodicImplicitTriangulation
    )

target_compile_definitions(triangulation
  PRIVATE
  TTK_IMPLICIT_PRECONDITIONS_THRESHOLD=${TTK_IMPLICIT_PRECONDITIONS_THRESHOLD}
)
