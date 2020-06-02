ttk_add_base_library(harmonicField
  SOURCES
    HarmonicField.cpp
  HEADERS
    HarmonicField.h
  DEPENDS
    geometry
    laplacian
    triangulation
  OPTIONAL_DEPENDS
    Eigen3::Eigen
    )

if (EIGEN3_FOUND)
  target_compile_definitions(harmonicField PUBLIC TTK_ENABLE_EIGEN)
endif()
