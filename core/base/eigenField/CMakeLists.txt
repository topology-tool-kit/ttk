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
  find_path(SPECTRA_INCLUDE_DIR Spectra/SymEigsSolver.h)
  if(SPECTRA_INCLUDE_DIR STREQUAL "SPECTRA_INCLUDE_DIR-NOTFOUND")
    option(TTK_ENABLE_SPECTRA "Enable Spectra support" OFF)
    message(STATUS "Spectra not found, disabling Spectra support in TTK.")
  else()
    option(TTK_ENABLE_SPECTRA "Enable Spectra support" ON)
  endif()
else()
  option(TTK_ENABLE_SPECTRA "Enable Spectra support" OFF)
  message(STATUS "Spectra not found, disabling Spectra support in TTK.")
endif()

if(TTK_ENABLE_SPECTRA)
  target_compile_definitions(eigenField PRIVATE TTK_ENABLE_SPECTRA)
  target_include_directories(eigenField SYSTEM PRIVATE ${SPECTRA_INCLUDE_DIR})
endif()
