ttk_add_base_library(eigenField
  SOURCES
    EigenField.cpp
  HEADERS
    EigenField.h
  LINK
    geometry
    laplacian
    triangulation
    )

if (EIGEN3_FOUND)
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
