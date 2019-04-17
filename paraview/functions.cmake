# Function to create the ParaView plugin for a TTK filter
# build one plugin per filter if TTK_BUILD_STANDALONE_PARAVIEW_PLUGINS is set

# call the _ttk_add_paraview_plugin function with the current folder as
# argument
macro(ttk_register_pv_filter pvFilter vtkModule)

  list(APPEND TTK_MODULES ${vtkModule})
  list(APPEND TTK_MODULE_FILES ${VTKWRAPPER_DIR}/${vtkModule}/vtk.module)

endmacro()
