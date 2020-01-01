# register a new filter to build in the TTK plugin
# deduce the location of the corresonding vtk.module file
# also register the xml file if given

set(DEBUG_WIDGETS "\
<IntVectorProperty name='Debug_UseAllCores' label='Use All Cores' command='SetUseAllCores' number_of_elements='1' default_values='1' panel_visibility='advanced'>\
    <BooleanDomain name='bool' />\
    <Documentation>Use all available cores.</Documentation>\
</IntVectorProperty>\
<IntVectorProperty name='Debug_ThreadNumber' label='Thread Number' command='SetThreadNumber' number_of_elements='1' default_values='1' panel_visibility='advanced'>\
    <IntRangeDomain name='range' min='1' max='256' />\
    <Hints>\
        <PropertyWidgetDecorator type='GenericDecorator' mode='visibility' property='Debug_UseAllCores' value='0' />\
    </Hints>\
    <Documentation>The maximum number of threads.</Documentation>\
</IntVectorProperty>\
<IntVectorProperty name='Debug_DebugLevel' label='Debug Level' command='SetDebugLevel' number_of_elements='1' default_values='3' panel_visibility='advanced'>\
    <IntRangeDomain name='range' min='0' max='5' />\
    <Documentation>Debug level.</Documentation>\
</IntVectorProperty>\
<Property name='Debug_Execute' label='Execute' command='Modified' panel_widget='command_button' panel_visibility='advanced'>\
    <Documentation>Executes the filter with the last applied parameters, which is handy to re-start pipeline execution from a specific element without changing parameters.</Documentation>\
</Property>\
<PropertyGroup panel_widget='Line' label='Testing'>\
    <Property name='Debug_UseAllCores' />\
    <Property name='Debug_ThreadNumber' />\
    <Property name='Debug_DebugLevel' />\
    <Property name='Debug_Execute' />\
</PropertyGroup>\
")

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
    if(NOT "${xmlFile}" STREQUAL "")

      # replace variables of original XML file and store generated file in build dir
      configure_file(${CMAKE_CURRENT_LIST_DIR}/${xmlFile} ${CMAKE_BINARY_DIR}/paraview/xmls/${xmlFile})

      # use generated xml file for pv plugin
      list(APPEND TTK_XMLS ${CMAKE_BINARY_DIR}/paraview/xmls/${xmlFile})
    endif()
  endif()
endmacro()
