ttk_add_base_library(harmonicField
  SOURCES
    HarmonicField.cpp
  HEADERS
    HarmonicField.h
  DEPENDS
    geometry
    laplacian
    triangulation
    )

if(TTK_ENABLE_EIGEN)
  target_compile_definitions(harmonicField PRIVATE TTK_ENABLE_EIGEN)
  target_link_libraries(harmonicField PRIVATE Eigen3::Eigen)
endif()
