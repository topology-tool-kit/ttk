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
    "${CMAKE_INSTALL_DATAROOTDIR}/ttk"
  )

target_compile_definitions(dimensionReduction PUBLIC -DTTK_SCRIPTS_PATH=${CMAKE_INSTALL_FULL_DATAROOTDIR}/ttk)

if(TTK_ENABLE_SCIKIT_LEARN)
  target_link_libraries(dimensionReduction PRIVATE Python3::Python Python3::NumPy)
endif()
