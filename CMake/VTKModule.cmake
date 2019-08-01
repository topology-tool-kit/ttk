# read the ttk.module file and create a list named _ttk_module_file from it
# this list can then be parsed using the cmake_parse_arguments method
macro(ttk_parse_module_file)
  # reconfigure when the module file is changed
  set_property(
    DIRECTORY
      "${CMAKE_CURRENT_SOURCE_DIR}"
    APPEND PROPERTY
    CMAKE_CONFIGURE_DEPENDS
      ${CMAKE_CURRENT_LIST_DIR}/ttk.module
  )
  # add_custom_target(CHECK_FILE_LISTS ALL ${CMAKE_CURRENT_LIST_DIR}/ttk.module)
  file(READ ${CMAKE_CURRENT_LIST_DIR}/ttk.module _ttk_module_file)
  # Replace comments.
  string(REGEX REPLACE "#[^\n]*\n" "\n" _ttk_module_file "${_ttk_module_file}")
  # Use argument splitting.
  string(REGEX REPLACE "( |\n)+" ";" _ttk_module_file "${_ttk_module_file}")
endmacro()

# add a new vtk module for the given name.
# also used to add the xml for paraview
macro(ttk_add_vtk_module)
  ttk_parse_module_file()
  cmake_parse_arguments("ARG" "" "NAME" "SOURCES;HEADERS;DEPENDS" ${_ttk_module_file})

  if(NOT TARGET ${ARG_NAME})
    vtk_module_add_module(${ARG_NAME}
      SOURCES
        ${ARG_SOURCES}
      HEADERS
        ${ARG_HEADERS}
      )
  endif()

  target_link_libraries(${ARG_NAME}
    PUBLIC
      ${VTK_LIBRARIES}
      ${ARG_DEPENDS}
    )
endmacro()

