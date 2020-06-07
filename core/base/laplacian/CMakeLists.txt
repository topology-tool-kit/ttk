ttk_add_base_library(laplacian
  SOURCES
    Laplacian.cpp
  HEADERS
    Laplacian.h
  DEPENDS
    geometry
    triangulation
  OPTIONAL_DEPENDS
    Eigen3::Eigen
    )

if (EIGEN3_FOUND)
  target_compile_definitions(laplacian PUBLIC TTK_ENABLE_EIGEN)
endif()
