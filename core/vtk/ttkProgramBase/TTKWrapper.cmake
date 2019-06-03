ttk_register_vtk_filter()

# not done by a cmake file here
target_link_libraries(ttkProgramBase PUBLIC ttkTriangulation)

if(MSVC)
  target_compile_definitions(ttkTriangulation PUBLIC vtkCommonDataModel_EXPORTS)
endif()
