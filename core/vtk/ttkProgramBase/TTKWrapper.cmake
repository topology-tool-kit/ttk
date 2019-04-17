ttk_register_vtk_filter()

if(MSVC)
  target_compile_definitions(ttkTriangulation PUBLIC vtkCommonDataModel_EXPORTS)
endif()
