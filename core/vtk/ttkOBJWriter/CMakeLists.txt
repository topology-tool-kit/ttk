ttk_add_vtk_module()

# TODO remove ?
if(MSVC)
  target_compile_definitions(ttkOBJWriter PUBLIC vtkIOLegacy_EXPORTS)
endif()
