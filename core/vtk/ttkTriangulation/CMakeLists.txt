ttk_add_vtk_library(ttkTriangulation
  SOURCES
    ttkTriangulation.cpp ttkWrapper.cpp
  HEADERS
    ttkTriangulation.h
    ttkWrapper.h
  LINK
    triangulation
    )

if (MSVC)
  target_compile_definitions(ttkTriangulation PUBLIC vtkCommonDataModel_EXPORTS)
endif()
