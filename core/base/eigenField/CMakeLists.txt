ttk_add_base_library(eigenField
  SOURCES
    EigenField.cpp
  HEADERS
    EigenField.h
  DEPENDS
    geometry
    laplacian
    triangulation
  OPTIONAL_DEPENDS
    Eigen3::Eigen
    )

if (EIGEN3_FOUND)
  target_compile_definitions(eigenField PUBLIC TTK_ENABLE_EIGEN)
endif()

if(TTK_ENABLE_SPECTRA)
  target_compile_definitions(eigenField PRIVATE TTK_ENABLE_SPECTRA)
  target_include_directories(eigenField SYSTEM PRIVATE ${SPECTRA_INCLUDE_DIR})
endif()
