# register a new filter to build in the TTK plugin
# deduce the location of the corresonding vtk.module file
# also register the xml file if given

file(READ "CMake/debug_widgets.xml" DEBUG_WIDGETS)
file(READ "CMake/topological_compression.xml" TOPOLOGICAL_COMPRESSION_WIDGETS)

macro(ttk_register_pv_filter vtkModuleDir xmlFile)
  if(NOT EXISTS "${VTKWRAPPER_DIR}/${vtkModuleDir}/vtk.module")
    message(FATAL_ERROR
      "Register a paraview module without the corresponding vtk.module: "
      ${VTKWRAPPER_DIR}/${vtkModuleDir}
      )
  endif()

  ttk_parse_module_file("${VTKWRAPPER_DIR}/${vtkModuleDir}/ttk.module")
  cmake_parse_arguments("TTK" "" "NAME" "SOURCES;HEADERS;DEPENDS" ${moduleFileContent})

  ttk_get_target(${TTK_NAME} TTK_TARGET)
  if(NOT "${TTK_TARGET}" STREQUAL "")
    list(APPEND TTK_MODULES ${TTK_TARGET})
    list(APPEND TTK_VTK_MODULE_FILES ${VTKWRAPPER_DIR}/${vtkModuleDir}/vtk.module)
    if(NOT "${xmlFile}" STREQUAL "")

      # replace variables of original XML file and store generated file in build dir
      configure_file(${CMAKE_CURRENT_LIST_DIR}/${xmlFile} ${CMAKE_BINARY_DIR}/paraview/xmls/${xmlFile})

      # use generated xml file for pv plugin
      list(APPEND TTK_XMLS ${CMAKE_BINARY_DIR}/paraview/xmls/${xmlFile})
    endif()
  endif()
endmacro()
