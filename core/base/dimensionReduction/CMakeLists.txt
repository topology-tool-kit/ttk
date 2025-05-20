ttk_add_base_library(dimensionReduction
  SOURCES
    DimensionReduction.cpp
  HEADERS
    DimensionReduction.h
  DEPENDS
    triangulation
    topoMap
    topologicallyConstrainedDimensionReduction
  )

install(
  FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/dimensionReduction.py
  DESTINATION
    ${TTK_SCRIPTS_PATH}
  )

# CPack doesn't accept absolute paths in install DESTINATION, must
# pass a relative path to TTK_SCRIPTS_PATH and append
# ${CMAKE_INSTALL_PREFIX} before giving it to dimensionReduction
if(IS_ABSOLUTE ${TTK_SCRIPTS_PATH})
  set(TTK_SCRIPTS_ABS_PATH "${TTK_SCRIPTS_PATH}")
else()
  set(TTK_SCRIPTS_ABS_PATH "${CMAKE_INSTALL_PREFIX}/${TTK_SCRIPTS_PATH}")
endif()

target_compile_definitions(dimensionReduction PUBLIC TTK_SCRIPTS_PATH=${TTK_SCRIPTS_ABS_PATH})

if(TTK_ENABLE_SCIKIT_LEARN)
  target_compile_definitions(dimensionReduction PUBLIC TTK_ENABLE_SCIKIT_LEARN)
  target_link_libraries(dimensionReduction PRIVATE Python3::Python Python3::NumPy)
endif()
