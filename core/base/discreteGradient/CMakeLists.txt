ttk_add_base_library(discreteGradient
  SOURCES
    DiscreteGradient.cpp
  HEADERS
    DiscreteGradient.h
    DiscreteGradient_Template.h
  LINK
    triangulation
    geometry
    scalarFieldCriticalPoints
    ftmTree
    )

option(TTK_ENABLE_DCG_OPTIMIZE_MEMORY "Enable Discrete Gradient memory optimization" OFF)
mark_as_advanced(TTK_ENABLE_DCG_OPTIMIZE_MEMORY)

if (TTK_ENABLE_DCG_OPTIMIZE_MEMORY)
  target_compile_definitions(discreteGradient PUBLIC TTK_ENABLE_DCG_OPTIMIZE_MEMORY)
endif()

