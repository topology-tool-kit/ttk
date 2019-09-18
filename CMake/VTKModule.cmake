# read the ttk.module file and create a list named moduleFileContent from it
# this list can then be parsed using the cmake_parse_arguments method
macro(ttk_parse_module_file moduleFile)
  # reconfigure when the module file is changed
  set_property(
    DIRECTORY
      "${CMAKE_CURRENT_SOURCE_DIR}"
    APPEND PROPERTY
    CMAKE_CONFIGURE_DEPENDS
      ${moduleFile}
  )
  # add_custom_target(CHECK_FILE_LISTS ALL ${CMAKE_CURRENT_LIST_DIR}/ttk.module)
  file(READ ${moduleFile} moduleFileContent)
  # Replace comments.
  string(REGEX REPLACE "#[^\n]*\n" "\n" moduleFileContent "${moduleFileContent}")
  # Use argument splitting.
  string(REGEX REPLACE "( |\n)+" ";" moduleFileContent "${moduleFileContent}")
endmacro()

# add a new vtk module for the given name.
# also used to add the xml for paraview
macro(ttk_add_vtk_module)
  ttk_parse_module_file(${CMAKE_CURRENT_LIST_DIR}/ttk.module)
  cmake_parse_arguments("ARG" "" "NAME" "SOURCES;HEADERS;DEPENDS" ${moduleFileContent})

  if(NOT TARGET ${ARG_NAME})
    vtk_module_add_module(${ARG_NAME}
      SOURCES
        ${ARG_SOURCES}
      HEADERS
        ${ARG_HEADERS}
      )
  endif()

  vtk_module_link(${ARG_NAME}
    PUBLIC
      ${VTK_LIBRARIES}
      ${ARG_DEPENDS}
    )
endmacro()

