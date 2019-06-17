# Function to create the library and install targets
# for a TTK BaseCode library.
#
# Usage:
# ttk_add_base_library(<library_name>
#     SOURCES <source list>
#     HEADERS <headers to install>
#     LINK <libraries to link>)
#
function(ttk_add_base_library library)
  cmake_parse_arguments(ARG "" "" "SOURCES;HEADERS;LINK" ${ARGN})

  set(CMAKE_SKIP_BUILD_RPATH TRUE)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/ttk/")
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

  add_library(${library} STATIC ${ARG_SOURCES} ${ARG_HEADERS})
  set_target_properties(${library} PROPERTIES
    POSITION_INDEPENDENT_CODE TRUE)

  target_include_directories(${library} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include/ttk/base>)

  if (ARG_LINK)
    target_link_libraries(${library} PUBLIC ${ARG_LINK})
  endif()

  ttk_set_compile_options(${library})

  install(TARGETS ${library}
    EXPORT TTKBaseTargets
    RUNTIME DESTINATION bin/ttk
    ARCHIVE DESTINATION lib/ttk
    LIBRARY DESTINATION lib/ttk)

  if (ARG_HEADERS)
    install(FILES ${ARG_HEADERS} DESTINATION include/ttk/base)
  endif()
endfunction()

# Function to create a header only template library and install targets
# for a TTK BaseCode template library. Note that here "libraries to link"
# has a bit of a different meaning. Since it's a template library no code
# is compiled or linked, but this "linking" is used to express dependencies
# on other libraries that users of this library will also need.
#
# Usage:
# ttk_add_base_template_library(<library_name>
#     SOURCES <headers list>
#     LINK <libraries to link>)
#
function(ttk_add_base_template_library library)
  cmake_parse_arguments(ARG "" "" "SOURCES;LINK" ${ARGN})

  add_library(${library} INTERFACE)
  target_include_directories(${library} INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include/ttk/base>)

  if(Boost_FOUND)
    target_link_libraries(${library} INTERFACE Boost::boost)
  endif()

  if (TTK_ENABLE_KAMIKAZE)
    target_compile_definitions(${library} INTERFACE TTK_ENABLE_KAMIKAZE)
  endif()

  string(TOUPPER "${CMAKE_BUILD_TYPE}" uppercase_CMAKE_BUILD_TYPE)
  if (uppercase_CMAKE_BUILD_TYPE MATCHES RELEASE)
    if (TTK_ENABLE_CPU_OPTIMIZATION)
      target_compile_options(${library} INTERFACE -march=native -O3)
    endif()
  endif()

  if (TTK_ENABLE_OPENMP)
    target_compile_definitions(${library} INTERFACE TTK_ENABLE_OPENMP)
    target_compile_options(${library} INTERFACE ${OpenMP_CXX_FLAGS})
  endif()

  if (ARG_LINK)
    target_link_libraries(${library} INTERFACE ${ARG_LINK})
  endif()

  if (TTK_ENABLE_MPI)
    target_link_libraries(${library} INTERFACE ${MPI_CXX_LIBRARIES})
  endif()

  if(TTK_ENABLE_EIGEN)
    target_link_libraries(${library} INTERFACE Eigen3::Eigen)
  endif()

  if (TTK_ENABLE_ZFP)
    target_link_libraries(${library} INTERFACE zfp::zfp)
  endif()

  if (TTK_ENABLE_GRAPHVIZ)
    target_compile_definitions(${library} INTERFACE TTK_ENABLE_GRAPHVIZ)
    target_include_directories(${library} INTERFACE ${GRAPHVIZ_INCLUDE_DIR})
  endif()

  if (TTK_ENABLE_SCIKIT_LEARN)
    target_compile_definitions(${library} INTERFACE TTK_ENABLE_SCIKIT_LEARN)
  endif()

  if (TTK_ENABLE_SQLITE3)
    target_compile_definitions(${library} INTERFACE TTK_ENABLE_SQLITE3)
    target_include_directories(${library} INTERFACE ${SQLITE3_INCLUDE_DIR})
  endif()

  if (GRAPHVIZ_INCLUDE_DIR)
    target_compile_definitions(${library} INTERFACE TTK_ENABLE_GRAPHVIZ)
    target_link_libraries(${library} INTERFACE ${GRAPHVIZ_CDT_LIBRARY})
    target_link_libraries(${library} INTERFACE ${GRAPHVIZ_GVC_LIBRARY})
    target_link_libraries(${library} INTERFACE ${GRAPHVIZ_CGRAPH_LIBRARY})
    target_link_libraries(${library} INTERFACE ${GRAPHVIZ_PATHPLAN_LIBRARY})
  endif()

  if (TTK_ENABLE_ZLIB)
    target_compile_definitions(${library} INTERFACE TTK_ENABLE_ZLIB)
    target_include_directories(${library} INTERFACE ${ZLIB_INCLUDE_DIR})
  endif()

  if (TTK_ENABLE_64BIT_IDS)
    target_compile_definitions(${library} INTERFACE TTK_ENABLE_64BIT_IDS)
  endif()

  install(TARGETS ${library}
    EXPORT TTKBaseTargets
    RUNTIME DESTINATION bin/ttk
    ARCHIVE DESTINATION lib/ttk
    LIBRARY DESTINATION lib/ttk)

  install(FILES ${ARG_SOURCES} DESTINATION include/ttk/base)
endfunction()
