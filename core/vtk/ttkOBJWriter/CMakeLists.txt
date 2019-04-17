ttk_add_vtk_library(ttkOBJWriter
  SOURCES
    ttkOBJWriter.cpp
  HEADERS
    ttkOBJWriter.h
    )

if (MSVC)
  target_compile_definitions(ttkOBJWriter PUBLIC vtkIOLegacy_EXPORTS)
endif()
