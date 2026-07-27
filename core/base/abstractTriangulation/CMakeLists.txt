ttk_add_base_library(abstractTriangulation
  SOURCES
    AbstractTriangulation.cpp
  HEADERS
    AbstractTriangulation.h
  DEPENDS
    common
    geometry
    )

option(TTK_ENABLE_DCG_OPTIMIZE_MEMORY "Enable Discrete Gradient memory optimization" OFF)
mark_as_advanced(TTK_ENABLE_DCG_OPTIMIZE_MEMORY)

if (TTK_ENABLE_DCG_OPTIMIZE_MEMORY)
  target_compile_definitions(abstractTriangulation PUBLIC TTK_ENABLE_DCG_OPTIMIZE_MEMORY)
endif()
