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

# add a new vtk module for with given name.
# This macro assume a ttk.module file is available
# in the current context
# it must have the following syntax:
# NAME: module name
# SOURCES: Specify the list of source files for the module
# HEADERS: Specify the list of header files to install for the module
# DEPENDS: Specify the link dependencies of the module
macro(ttk_add_vtk_module)
  ttk_parse_module_file(${CMAKE_CURRENT_LIST_DIR}/ttk.module)
  cmake_parse_arguments("TTK" "" "NAME" "SOURCES;HEADERS;DEPENDS" ${moduleFileContent})

  if(NOT TARGET ${TTK_NAME})
    vtk_module_add_module(${TTK_NAME}
      SOURCES
        ${TTK_SOURCES}
      HEADERS
        ${TTK_HEADERS}
      )
  endif()

  vtk_module_link(${TTK_NAME}
    PUBLIC
      ${VTK_LIBRARIES}
      ${TTK_DEPENDS}
    )

  # Fix a race condition in the VTK's CMake:
  # https://discourse.vtk.org/t/building-vtk-modules-with-dependencies-results-in-race-condition-in-make/1711
  if(TARGET ${TTK_NAME}-hierarchy)
    add_dependencies(${TTK_NAME} ${TTK_NAME}-hierarchy)
  endif()
endmacro()

# deal with whitelist mechanism
# return the target name if target is enabled
# empty otherwise
function(ttk_get_target ttk_module ttk_target)
  string(TOUPPER ${ttk_module} TTK_UPPER_NAME)
  if(NOT DEFINED TTK_BUILD_${TTK_UPPER_NAME})
    option(TTK_BUILD_${TTK_UPPER_NAME} "Build the ${TTK_UPPER_NAME} filter" ${TTK_ENABLE_FILTER_DEFAULT})
    mark_as_advanced(TTK_BUILD_${TTK_UPPER_NAME})
  endif()
  if(${TTK_BUILD_${TTK_UPPER_NAME}})
    set(${ttk_target} ${ttk_module} PARENT_SCOPE)
  else()
    message(STATUS "Disabled: ${ttk_module}")
    set(${ttk_target} "" PARENT_SCOPE)
  endif()
endfunction()
