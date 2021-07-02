# Function to create the library and install targets
# for a TTK BaseCode library.
#
# Usage:
# ttk_add_base_library(<library_name>
#     SOURCES <source list>
#     HEADERS <headers to install>
#     DEPENDS <libraries to link>)
#
function(ttk_add_base_library library)
  cmake_parse_arguments(ARG "" "" "SOURCES;HEADERS;DEPENDS;OPTIONAL_DEPENDS" ${ARGN})

  if(TTK_ENABLE_SHARED_BASE_LIBRARIES)
    add_library(${library} SHARED ${ARG_SOURCES} ${ARG_HEADERS})
  else()
    add_library(${library} STATIC ${ARG_SOURCES} ${ARG_HEADERS})
  endif()

  set_target_properties(${library} PROPERTIES OUTPUT_NAME "ttkBase${library}")

  set_target_properties(${library}
    PROPERTIES
      POSITION_INDEPENDENT_CODE TRUE
    )
  target_include_directories(${library}
    PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/ttk/base>
    )

  if(ARG_DEPENDS)
    target_link_libraries(${library} PUBLIC ${ARG_DEPENDS})
  endif()

  if(ARG_OPTIONAL_DEPENDS)
    foreach(OPT_DEP IN LISTS ARG_OPTIONAL_DEPENDS)
      if(TARGET ${OPT_DEP})
       target_link_libraries(${library} PUBLIC ${OPT_DEP})
      endif()
    endforeach()
  endif()

  ttk_set_compile_options(${library})

  install(TARGETS ${library}
    EXPORT TTKBaseTargets
    )

  if(ARG_HEADERS)
    install(FILES ${ARG_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ttk/base)
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
#     HEADERS <headers list>
#     DEPENDS <libraries to link>)
#
function(ttk_add_base_template_library library)
  cmake_parse_arguments(ARG "" "" "HEADERS;DEPENDS" ${ARGN})

  add_library(${library} INTERFACE)

  target_include_directories(${library}
    INTERFACE
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/ttk/base>
    )

  if(TTK_ENABLE_DOUBLE_TEMPLATING)
    target_compile_definitions(${library} INTERFACE TTK_ENABLE_DOUBLE_TEMPLATING)
  endif()

  if(TTK_ENABLE_KAMIKAZE)
    target_compile_definitions(${library} INTERFACE TTK_ENABLE_KAMIKAZE)
  endif()

  string(TOUPPER "${CMAKE_BUILD_TYPE}" uppercase_CMAKE_BUILD_TYPE)
  if(uppercase_CMAKE_BUILD_TYPE MATCHES RELEASE)
    if(TTK_ENABLE_CPU_OPTIMIZATION AND NOT MSVC)
      target_compile_options(${library} INTERFACE -march=native -O3)
    endif()
  endif()

  if(TTK_ENABLE_OPENMP)
    target_compile_definitions(${library} INTERFACE TTK_ENABLE_OPENMP)
    target_link_libraries(${library} INTERFACE OpenMP::OpenMP_CXX)
  endif()

  if(ARG_DEPENDS)
    target_link_libraries(${library} INTERFACE ${ARG_DEPENDS})
  endif()

  if(TTK_ENABLE_MPI)
    target_link_libraries(${library} INTERFACE ${MPI_CXX_LIBRARIES})
  endif()

  if(TTK_ENABLE_64BIT_IDS)
    target_compile_definitions(${library} INTERFACE TTK_ENABLE_64BIT_IDS)
  endif()

  install(TARGETS ${library}
    EXPORT TTKBaseTargets
    )

  install(FILES ${ARG_HEADERS} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/ttk/base")
endfunction()

# Add compile flags and defintions to the target
# according to the options selected by the user.
# Only options and defintions common to all modules
# should be here.
#
# Usage:
# ttk_set_compile_options(<library_name>)
#
function(ttk_set_compile_options library)

  # compilation flags
  target_compile_options(${library} PRIVATE ${TTK_COMPILER_FLAGS})

  # linker flags
  if(TTK_LINKER_FLAGS)
    target_link_options(${library} PRIVATE ${TTK_LINKER_FLAGS})
  endif()

  if (TTK_ENABLE_DOUBLE_TEMPLATING)
    target_compile_definitions(${library} PUBLIC TTK_ENABLE_DOUBLE_TEMPLATING)
  endif()

  if (TTK_ENABLE_SHARED_BASE_LIBRARIES)
    target_compile_definitions(${library} PUBLIC TTK_ENABLE_SHARED_BASE_LIBRARIES)
  endif()

  if (TTK_ENABLE_KAMIKAZE)
    target_compile_definitions(${library} PUBLIC TTK_ENABLE_KAMIKAZE)
  endif()

  if (TTK_ENABLE_OPENMP)
    target_compile_definitions(${library} PUBLIC TTK_ENABLE_OPENMP)
    target_link_libraries(${library} PUBLIC OpenMP::OpenMP_CXX)
  endif()

  if (TTK_ENABLE_MPI)
    target_compile_definitions(${library} PUBLIC TTK_ENABLE_MPI)
    target_include_directories(${library} PUBLIC ${MPI_CXX_INCLUDE_PATH})
    target_link_libraries(${library} PUBLIC ${MPI_CXX_LIBRARIES})
  endif()

  if (TTK_ENABLE_64BIT_IDS)
    target_compile_definitions(${library} PUBLIC TTK_ENABLE_64BIT_IDS)
  endif()

endfunction()
