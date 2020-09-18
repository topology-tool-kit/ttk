ttk_add_base_library(harmonicField
  SOURCES
    HarmonicField.cpp
  HEADERS
    HarmonicField.h
  DEPENDS
    geometry
    laplacian
    triangulation
    )

if (EIGEN3_FOUND)
  target_compile_definitions(harmonicField PRIVATE TTK_ENABLE_EIGEN)
  target_include_directories(harmonicField SYSTEM PRIVATE ${EIGEN3_INCLUDE_DIR})
endif()
