# register a new filter to build in the TTK plugin
# also deduce the location of the corresonding vtk.module file
macro(ttk_register_pv_filter pvFilter vtkModule)
  if(NOT EXISTS ${VTKWRAPPER_DIR}/${vtkModule}/vtk.module)
    message(FATAL_ERROR
      "Register a paraview filter without the corresponding vtk.module: "
      ${VTKWRAPPER_DIR}/${vtkModule}
      )
  endif()
  list(APPEND TTK_MODULES ${vtkModule})
  list(APPEND TTK_MODULE_FILES ${VTKWRAPPER_DIR}/${vtkModule}/vtk.module)
endmacro()
