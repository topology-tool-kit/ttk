ttk_add_vtk_module()

# TODO remove ?
if(MSVC)
  target_compile_definitions(ttkOFFWriter PUBLIC vtkIOLegacy_EXPORTS)
endif()
