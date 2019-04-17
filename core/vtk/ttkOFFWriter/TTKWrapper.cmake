ttk_register_vtk_filter()

if(MSVC)
  target_compile_definitions(ttkOFFWriter PUBLIC vtkIOLegacy_EXPORTS)
endif()
