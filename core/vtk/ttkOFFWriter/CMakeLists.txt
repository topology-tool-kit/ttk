ttk_add_vtk_library(ttkOFFWriter
  SOURCES
    ttkOFFWriter.cpp
  HEADERS
    ttkOFFWriter.h
    )

if (MSVC)
  target_compile_definitions(ttkOFFWriter PUBLIC vtkIOLegacy_EXPORTS)
endif()
