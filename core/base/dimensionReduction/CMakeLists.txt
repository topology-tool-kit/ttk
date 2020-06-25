ttk_find_python()

ttk_add_base_library(dimensionReduction
  SOURCES
    DimensionReduction.cpp
  HEADERS
    DimensionReduction.h
  DEPENDS
    triangulation
  )

install(
  FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/dimensionReduction.py
  DESTINATION
    ${CMAKE_INSTALL_PREFIX}/share/scripts/ttk
  )

target_compile_definitions(dimensionReduction PUBLIC -DTTK_SCRIPTS_PATH=${CMAKE_INSTALL_PREFIX}/share/scripts/ttk)
