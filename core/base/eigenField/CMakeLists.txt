ttk_add_base_library(eigenField
  SOURCES
    EigenField.cpp
  HEADERS
    EigenField.h
  DEPENDS
    geometry
    laplacian
    triangulation
    )

if(TTK_ENABLE_SPECTRA)
  target_compile_definitions(eigenField PRIVATE TTK_ENABLE_EIGEN)
  target_include_directories(eigenField SYSTEM PRIVATE ${EIGEN3_INCLUDE_DIR})
  target_compile_definitions(eigenField PRIVATE TTK_ENABLE_SPECTRA)
  target_link_libraries(eigenField PRIVATE Spectra::Spectra)
endif()
