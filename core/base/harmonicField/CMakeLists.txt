ttk_add_base_library(harmonicField
  SOURCES
    HarmonicField.cpp
  HEADERS
    HarmonicField.h
  LINK
    geometry
    laplacian
    triangulation
    )

target_include_directories(harmonicField PRIVATE "/home/pierre/Projects/extern/spectra/include")
#target_compile_options(harmonicField PRIVATE "-O1")
