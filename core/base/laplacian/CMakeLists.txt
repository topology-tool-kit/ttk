ttk_add_base_library(laplacian
  SOURCES
    Laplacian.cpp
  HEADERS
    Laplacian.h
  DEPENDS
    geometry
    triangulation
    )

if(TTK_ENABLE_EIGEN)
  target_compile_definitions(laplacian PRIVATE TTK_ENABLE_EIGEN)
  target_link_libraries(laplacian PRIVATE Eigen3::Eigen)
endif()
