ttk_register_vtk_filter()

if(MSVC)
  target_compile_definitions(ttkOBJWriter PUBLIC vtkIOLegacy_EXPORTS)
endif()
