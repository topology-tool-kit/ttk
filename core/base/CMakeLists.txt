include(functions.cmake)

# traverse subdirectories
set(LIB_LIST "")
file(GLOB BASECODE_DIRS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *)
foreach(BASECODE ${BASECODE_DIRS})
  if (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${BASECODE})
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/${BASECODE})
    set(LIB_LIST "${BASECODE};${LIB_LIST}")
  endif()
endforeach()

add_library(baseAll INTERFACE)
target_link_libraries(baseAll INTERFACE ${LIB_LIST})

install(TARGETS baseAll
  EXPORT TTKBaseTargets
  RUNTIME DESTINATION bin/ttk
  ARCHIVE DESTINATION lib/ttk
  LIBRARY DESTINATION lib/ttk)

install(EXPORT TTKBaseTargets
  NAMESPACE ttk::base::
  DESTINATION lib/cmake/ttk)

configure_file(TTKBaseConfig.cmake.in TTKBaseConfig.cmake @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/TTKBaseConfig.cmake"
  DESTINATION lib/cmake/ttk)

