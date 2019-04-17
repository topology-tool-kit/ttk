# Function to create the library and install targets
# for a TTK VTK Wrapper library. The VTK Libraries will
# automatically be linked.
# Options:
# SOURCES: Specify the list of source files for the library
# HEADERS: Specify the list of header files to install for the library
# DEPENDS: Specify the link dependencies of the library
#
# To build a VTK Wrapper library:
# ttk_add_vtk_library(<library_name>
#     SOURCES <source list>
#     HEADERS <headers to install>
#     DEPENDS <libraries to link>)
#
function(_ttk_add_vtk_library library)
  cmake_parse_arguments(ARG "" "CUR_FOLD" "SOURCES;HEADERS;DEPENDS" ${ARGN})

  set(CMAKE_SKIP_BUILD_RPATH TRUE)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/ttk/")
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

  # Build the VTK Wrapper
  list(TRANSFORM ARG_SOURCES PREPEND "${ARG_CUR_FOLD}/")
  add_library(${library} SHARED ${ARG_SOURCES})
  target_link_libraries(${library} PUBLIC ${VTK_LIBRARIES} ${ARG_DEPENDS})

  target_include_directories(${library} PUBLIC
    $<BUILD_INTERFACE:${ARG_CUR_FOLD}>
    $<INSTALL_INTERFACE:include/ttk/vtk>)

  ttk_set_compile_options(${library})

  if(MSVC)
    target_compile_definitions(${library} PUBLIC vtkFiltersCore_EXPORTS)
  endif()

  install(TARGETS ${library}
    EXPORT TTKVTKTargets
    RUNTIME DESTINATION bin/ttk
    ARCHIVE DESTINATION lib/ttk
    LIBRARY DESTINATION lib/ttk)

  if(ARG_HEADERS)
    list(TRANSFORM ARG_HEADERS PREPEND "${ARG_CUR_FOLD}/")
    install(FILES ${ARG_HEADERS} DESTINATION include/ttk/vtk)
  endif()
endfunction()

# TODO remove me
macro(ttk_add_vtk_library library)
  _ttk_add_vtk_library(${library}
    ${ARGN}
    CUR_FOLD
      ${CMAKE_CURRENT_LIST_DIR}
      )
endmacro()

# read the ttk.module file and create a list named _ttk_module_file from it
macro(ttk_parse_module_file)
  # reconfigure when this file is changed
  set_property(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}" APPEND
    PROPERTY
      CMAKE_CONFIGURE_DEPENDS ${CMAKE_CURRENT_LIST_DIR}/ttk.module)
  # add_custom_target(CHECK_FILE_LISTS ALL ${CMAKE_CURRENT_LIST_DIR}/ttk.module)
  file(READ ${CMAKE_CURRENT_LIST_DIR}/ttk.module _ttk_module_file)
  # Replace comments.
  string(REGEX REPLACE "#[^\n]*\n" "\n" _ttk_module_file "${_ttk_module_file}")
  # Use argument splitting.
  string(REGEX REPLACE "( |\n)+" ";" _ttk_module_file "${_ttk_module_file}")
endmacro()

# Use the ttk.module file and call the ttk_add_library function with the
# current dir
macro(ttk_register_vtk_filter)
  ttk_parse_module_file()
  cmake_parse_arguments( "ARG" "" "NAME" "SOURCES;HEADERS;DEPENDS;XMLS" ${_ttk_module_file})

  _ttk_add_vtk_library(${ARG_NAME}
  SOURCES
    ${ARG_SOURCES}
  HEADERS
    ${ARG_HEADERS}
  DEPENDS
    ${ARG_DEPENDS}
  CUR_FOLD
    ${CMAKE_CURRENT_LIST_DIR}
    )
endmacro()

macro(ttk_add_vtk_module)
  ttk_parse_module_file()
  cmake_parse_arguments( "ARG" "" "NAME" "SOURCES;HEADERS;DEPENDS;XMLS" ${_ttk_module_file})

  vtk_module_add_module(${ARG_NAME}
    SOURCES
      ${ARG_SOURCES}
    HEADERS
      ${ARG_HEADERS}
      )

  target_link_libraries(${ARG_NAME}
    PUBLIC
      ${ARG_DEPENDS}
      )

  paraview_add_server_manager_xmls(
    XMLS
      ${ARG_XMLS}
    )
endmacro()


