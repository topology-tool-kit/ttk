
# When building ParaView plugins, we need to
# handle installation of header files by hand
macro(ttk_module_register_headers filterPath)
  ttk_parse_module_file(${filterPath}/ttk.module)
  cmake_parse_arguments("TTK" "" "NAME" "SOURCES;HEADERS;DEPENDS;XMLS" ${moduleFileContent})

  if(TTK_HEADERS)
    vtk_module_install_headers(
      FILES
        ${TTK_HEADERS}
      SUBDIR
        include/ttk
      )
  endif()

endmacro()
