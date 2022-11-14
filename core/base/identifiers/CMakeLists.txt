ttk_add_base_library(identifiers
  SOURCES
    Identifiers.cpp
  HEADERS
    Identifiers.h
  DEPENDS
    triangulation
    kdTree
)
target_compile_options(identifiers PUBLIC -O0 -g)
