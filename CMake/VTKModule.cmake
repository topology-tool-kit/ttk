# read the ttk.module file and create a list named _ttkModuleFileContent from it
# this list can then be parsed using the cmake_parse_arguments method
macro(ttk_parse_module_file ttkModuleFile)
  # reconfigure when the module file is changed
  set_property(
    DIRECTORY
      "${CMAKE_CURRENT_SOURCE_DIR}"
    APPEND PROPERTY
    CMAKE_CONFIGURE_DEPENDS
      ${ttkModuleFile}
  )
  # add_custom_target(CHECK_FILE_LISTS ALL ${CMAKE_CURRENT_LIST_DIR}/ttk.module)
  file(READ ${ttkModuleFile} _ttkModuleFileContent)
  # Replace comments.
  string(REGEX REPLACE "#[^\n]*\n" "\n" _ttkModuleFileContent "${_ttkModuleFileContent}")
  # Use argument splitting.
  string(REGEX REPLACE "( |\n)+" ";" _ttkModuleFileContent "${_ttkModuleFileContent}")
endmacro()

# add a new vtk module for the given name.
# also used to add the xml for paraview
macro(ttk_add_vtk_module)
  ttk_parse_module_file(${CMAKE_CURRENT_LIST_DIR}/ttk.module)
  cmake_parse_arguments("ARG" "" "NAME" "SOURCES;HEADERS;DEPENDS" ${_ttkModuleFileContent})

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

