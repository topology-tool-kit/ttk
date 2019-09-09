ttk_find_python()

ttk_add_base_library(dimensionReduction
  SOURCES
    DimensionReduction.cpp
  HEADERS
    DimensionReduction.h
  LINK
    triangulation
    )

install(
  FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/dimensionReduction.py
  DESTINATION
    ${CMAKE_INSTALL_PREFIX}/share/scripts/ttk
    )

add_definitions(-DTTK_SCRIPTS_PATH=${CMAKE_INSTALL_PREFIX}/share/scripts/ttk)
