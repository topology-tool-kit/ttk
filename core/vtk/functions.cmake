# Function to create the library and install targets
# for a TTK VTK Wrapper library. The VTK Libraries will
# automatically be linked.
# Options:
# SOURCES: Specify the list of source files for the library
# HEADERS: Specify the list of header files to install for the library
# LINK: Specify the link dependencies of the library
#
# To build a VTK Wrapper library:
# ttk_add_vtk_library(<library_name>
#     SOURCES <source list>
#     HEADERS <headers to install>
#     LINK <libraries to link>)
#

function(ttk_add_vtk_library library)
  cmake_parse_arguments(ARG "" "" "SOURCES;HEADERS;LINK" ${ARGN})

  set(CMAKE_SKIP_BUILD_RPATH TRUE)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/ttk/")
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

  # Build the VTK Wrapper
  add_library(${library} SHARED ${ARG_SOURCES})
  target_link_libraries(${library} PUBLIC ${VTK_LIBRARIES} ${ARG_LINK})
  # add HEADERS paramter to target properties
  set_target_properties(${library} PROPERTIES TTK_HEADERS "${ARG_HEADERS}")

  target_include_directories(${library} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include/ttk/vtk>)

  ttk_set_compile_options(${library})

  if (MSVC)
    target_compile_definitions(${library} PUBLIC vtkFiltersCore_EXPORTS)
  endif()

  install(TARGETS ${library}
    EXPORT TTKVTKTargets
    RUNTIME DESTINATION bin/ttk
    ARCHIVE DESTINATION lib/ttk
    LIBRARY DESTINATION lib/ttk)

  if (ARG_HEADERS)
    install(FILES ${ARG_HEADERS} DESTINATION include/ttk/vtk)
  endif()
endfunction()

