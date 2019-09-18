# register a new filter to build in the TTK plugin
# also deduce the location of the corresonding vtk.module file
macro(ttk_register_pv_filter vtkModuleDir xmlFile)
  if(NOT EXISTS ${VTKWRAPPER_DIR}/${vtkModuleDir}/vtk.module)
    message(FATAL_ERROR
      "Register a paraview filter without the corresponding vtk.module: "
      ${VTKWRAPPER_DIR}/${vtkModuleDir}
      )
  endif()
  list(APPEND TTK_XMLS ${CMAKE_CURRENT_LIST_DIR}/${xmlFile})
  list(APPEND TTK_MODULES ${vtkModuleDir})
  list(APPEND TTK_VTK_MODULE_FILESS ${VTKWRAPPER_DIR}/${vtkModuleDir}/vtk.module)
endmacro()

# for filters without xml
macro(ttk_register_pv_module vtkModuleDir)
  if(NOT EXISTS ${VTKWRAPPER_DIR}/${vtkModuleDir}/vtk.module)
    message(FATAL_ERROR
      "Register a paraview module without the corresponding vtk.module: "
      ${VTKWRAPPER_DIR}/${vtkModuleDir}
      )
  endif()
  list(APPEND TTK_MODULES ${vtkModuleDir})
  list(APPEND TTK_VTK_MODULE_FILESS ${VTKWRAPPER_DIR}/${vtkModuleDir}/vtk.module)
endmacro()
