# for the ttk_parse_module_file
include(CMake/VTKModule.cmake)

# register a new filter to build in the TTK plugin
# deduce the location of the corresponding vtk.module file
# also register the xml file if given

# TODO ... this has nothing to do here.
file(READ "CMake/debug_widgets.xml" DEBUG_WIDGETS)
file(READ "CMake/topological_compression.xml" TOPOLOGICAL_COMPRESSION_WIDGETS)
file(READ "CMake/merge_tree_input.xml" MERGE_TREE_INPUT_WIDGETS)
file(READ "CMake/merge_tree_preprocess.xml" MERGE_TREE_PREPROCESS_WIDGETS)
file(READ "CMake/merge_tree_planar_layout.xml" MERGE_TREE_PLANAR_LAYOUT_WIDGETS)

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

function(ttk_set_paraview_rpath TARGET_NAME)
  if(APPLE)
    # On macOS,
    # look into the subdirectory "TopologyToolKit"
    # to find the actual plugins.
    get_target_property(TEMP
        ${TARGET_NAME} INSTALL_RPATH
        )
    set_target_properties(${TARGET_NAME}
      PROPERTIES
        INSTALL_RPATH "@loader_path/${TTK_PLUGIN_SUBDIR};${TEMP}"
    )
  endif(APPLE)
endfunction(ttk_set_paraview_rpath TARGET_NAME)

