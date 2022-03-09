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
    ${TTK_SCRIPTS_PATH}
  )

target_compile_definitions(dimensionReduction PUBLIC TTK_SCRIPTS_PATH=${TTK_SCRIPTS_PATH})

if(TTK_ENABLE_SCIKIT_LEARN)
  target_compile_definitions(dimensionReduction PRIVATE TTK_ENABLE_SCIKIT_LEARN)
  target_link_libraries(dimensionReduction PRIVATE Python3::Python Python3::NumPy)
endif()
