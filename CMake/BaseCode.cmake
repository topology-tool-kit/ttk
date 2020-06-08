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

  add_library(${library}
    STATIC
      ${ARG_SOURCES}
      ${ARG_HEADERS}
    )
  set_target_properties(${library}
    PROPERTIES
      POSITION_INDEPENDENT_CODE TRUE
    )
  target_include_directories(${library}
    PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
      $<INSTALL_INTERFACE:include/ttk/base>
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
    EXPORT
      TTKBaseTargets
    RUNTIME DESTINATION
      "${CMAKE_INSTALL_BINDIR}/ttk"
    ARCHIVE DESTINATION
      "${CMAKE_INSTALL_LIBDIR}/ttk"
    LIBRARY DESTINATION
      "${CMAKE_INSTALL_LIBDIR}/ttk"
    )

  if(ARG_HEADERS)
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
#     HEADERS <headers list>
#     DEPENDS <libraries to link>)
#
function(ttk_add_base_template_library library)
  cmake_parse_arguments(ARG "" "" "HEADERS;DEPENDS" ${ARGN})

  add_library(${library} INTERFACE)

  target_include_directories(${library}
    INTERFACE
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
      $<INSTALL_INTERFACE:include/ttk/base>
    )

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
    target_compile_options(${library} INTERFACE ${OpenMP_CXX_FLAGS})
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
    RUNTIME DESTINATION
      "${CMAKE_INSTALL_BINDIR}/ttk"
    ARCHIVE DESTINATION
      "${CMAKE_INSTALL_LIBDIR}/ttk"
    LIBRARY DESTINATION
      "${CMAKE_INSTALL_LIBDIR}/ttk"
    )

  install(FILES ${ARG_HEADERS} DESTINATION include/ttk/base)
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

  if (TTK_ENABLE_KAMIKAZE)
    target_compile_definitions(${library} PUBLIC TTK_ENABLE_KAMIKAZE)
  endif()

  if (TTK_ENABLE_OPENMP)
    target_compile_definitions(${library} PUBLIC TTK_ENABLE_OPENMP)
    target_compile_options(${library} PUBLIC ${OpenMP_CXX_FLAGS})
    target_link_libraries(${library} PUBLIC ${OpenMP_CXX_LIBRARIES})

    if (TTK_ENABLE_OMP_PRIORITY)
      target_compile_definitions(${library} PUBLIC TTK_ENABLE_OMP_PRIORITY)
    endif()
  endif()

  if (TTK_ENABLE_MPI)
    target_compile_definitions(${library} PUBLIC TTK_ENABLE_MPI)
    target_include_directories(${library} PUBLIC ${MPI_CXX_INCLUDE_PATH})
    target_link_libraries(${library} PUBLIC ${MPI_CXX_LIBRARIES})
  endif()

  if (TTK_ENABLE_SCIKIT_LEARN)
    target_compile_definitions(${library} PUBLIC TTK_ENABLE_SCIKIT_LEARN)
  endif()

  if (TTK_ENABLE_ZLIB)
    target_compile_definitions(${library} PUBLIC TTK_ENABLE_ZLIB)
  endif()

  # TODO per module
  if (TTK_ENABLE_GRAPHVIZ AND GRAPHVIZ_FOUND)
    target_compile_definitions(${library} PUBLIC TTK_ENABLE_GRAPHVIZ)
    target_include_directories(${library} PUBLIC ${GRAPHVIZ_INCLUDE_DIR})
    target_link_libraries(${library} PUBLIC ${GRAPHVIZ_CDT_LIBRARY})
    target_link_libraries(${library} PUBLIC ${GRAPHVIZ_GVC_LIBRARY})
    target_link_libraries(${library} PUBLIC ${GRAPHVIZ_CGRAPH_LIBRARY})
    target_link_libraries(${library} PUBLIC ${GRAPHVIZ_PATHPLAN_LIBRARY})
  endif()

  # TODO per module
  if (TTK_ENABLE_SQLITE3)
    target_compile_definitions(${library} PUBLIC TTK_ENABLE_SQLITE3)
    target_include_directories(${library} PUBLIC ${SQLITE3_INCLUDE_DIR})
    target_link_libraries(${library} PUBLIC ${SQLITE3_LIBRARY})
  endif()

  if (TTK_ENABLE_64BIT_IDS)
    target_compile_definitions(${library} PUBLIC TTK_ENABLE_64BIT_IDS)
  endif()

endfunction()

# Used by basedCode requiring "Python.h"

function(ttk_find_python)
  find_package(PythonLibs QUIET)

  if(PYTHON_INCLUDE_DIRS)
    include_directories(SYSTEM ${PYTHON_INCLUDE_DIRS})

    if (${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.15")
      string(REPLACE \".\" \" \"
        PYTHON_VERSION_LIST ${PYTHONLIBS_VERSION_STRING})
    else()
      string(REPLACE "." " "
        PYTHON_VERSION_LIST ${PYTHONLIBS_VERSION_STRING})
    endif()
    separate_arguments(PYTHON_VERSION_LIST)
    list(GET PYTHON_VERSION_LIST 0 PYTHON_MAJOR_VERSION)
    list(GET PYTHON_VERSION_LIST 1 PYTHON_MINOR_VERSION)

    set(TTK_PYTHON_MAJOR_VERSION "${PYTHON_MAJOR_VERSION}"
      CACHE INTERNAL "TTK_PYTHON_MAJOR_VERSION")
    set(TTK_PYTHON_MINOR_VERSION "${PYTHON_MINOR_VERSION}"
      CACHE INTERNAL "TTK_PYTHON_MINOR_VERSION")

    if(TTK_PYTHON_MAJOR_VERSION)
      message(STATUS "Python version: ${TTK_PYTHON_MAJOR_VERSION}.${TTK_PYTHON_MINOR_VERSION}")
    else()
      message(STATUS "Python version: NOT-FOUND")
    endif()

    find_path(PYTHON_NUMPY_INCLUDE_DIR numpy/arrayobject.h PATHS
      ${PYTHON_INCLUDE_DIRS}
      /usr/lib/python${PYTHON_MAJOR_VERSION}.${PYTHON_MINOR_VERSION}/site-packages/numpy/core/include/
      /usr/local/lib/python${PYTHON_MAJOR_VERSION}.${PYTHON_MINOR_VERSION}/site-packages/numpy/core/include)
    if(PYTHON_NUMPY_INCLUDE_DIR)
      message(STATUS "Numpy headers: ${PYTHON_NUMPY_INCLUDE_DIR}")
      include_directories(SYSTEM ${PYTHON_NUMPY_INCLUDE_DIR})
    else()
      message(STATUS "Numpy headers: NOT-FOUND")
    endif()
  endif()

  if(PYTHON_NUMPY_INCLUDE_DIR)
    option(TTK_ENABLE_SCIKIT_LEARN "Enable scikit-learn support" ON)
  else()
    option(TTK_ENABLE_SCIKIT_LEARN "Enable scikit-learn support" OFF)
    message(STATUS 
      "Improper python/numpy setup. Disabling sckikit-learn support in TTK.")
  endif()

endfunction()
