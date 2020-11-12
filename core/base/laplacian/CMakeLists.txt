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
  target_include_directories(laplacian SYSTEM PRIVATE ${EIGEN3_INCLUDE_DIR})
endif()
