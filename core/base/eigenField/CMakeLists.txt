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

target_include_directories(eigenField PRIVATE "/home/pierre/Projects/extern/spectra/include")
