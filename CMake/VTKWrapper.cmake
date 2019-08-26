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

  if(NOT TARGET ${library})
    # Build the VTK Wrapper
    list(TRANSFORM ARG_SOURCES PREPEND "${ARG_CUR_FOLD}/")
    add_library(${library} SHARED ${ARG_SOURCES})

    target_link_libraries(${library}
      PUBLIC
        ${VTK_LIBRARIES}
        ${ARG_DEPENDS}
    )

    target_include_directories(${library}
      PUBLIC
        $<BUILD_INTERFACE:${ARG_CUR_FOLD}>
        $<INSTALL_INTERFACE:include/ttk/vtk>
    )

    ttk_set_compile_options(${library})
  endif()

  if(MSVC)
    target_compile_definitions(${library} PUBLIC vtkFiltersCore_EXPORTS)
  endif()

  include(GNUInstallDirs)
  install(
    TARGETS
      ${library}
    EXPORT
      TTKVTKTargets
    RUNTIME DESTINATION
      ${CMAKE_INSTALL_BINDIR}/ttk
    ARCHIVE DESTINATION
      ${CMAKE_INSTALL_LIBDIR}/ttk
    LIBRARY DESTINATION
      ${CMAKE_INSTALL_LIBDIR}/ttk
  )

  if(ARG_HEADERS)
    list(TRANSFORM ARG_HEADERS PREPEND "${ARG_CUR_FOLD}/")
    install(FILES ${ARG_HEADERS} DESTINATION include/ttk/vtk)
  endif()
endfunction()

# Use the ttk.module file and call the ttk_add_library function with the
# current dir
macro(ttk_register_vtk_filter filterPath)
  ttk_parse_module_file(${filterPath}/ttk.module)
  cmake_parse_arguments( "ARG" "" "NAME" "SOURCES;HEADERS;DEPENDS;XMLS" ${_ttkModuleFileContent})

  _ttk_add_vtk_library(${ARG_NAME}
    SOURCES
      ${ARG_SOURCES}
    HEADERS
      ${ARG_HEADERS}
    DEPENDS
      ${ARG_DEPENDS}
    CUR_FOLD
      ${filterPath}
  )
endmacro()
