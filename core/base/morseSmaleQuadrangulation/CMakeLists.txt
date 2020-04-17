if ("${TTK_CELL_ARRAY_LAYOUT}" STREQUAL "SingleArray")
  ttk_add_base_library(morseSmaleQuadrangulation
    SOURCES
      MorseSmaleQuadrangulation.cpp
    HEADERS
      MorseSmaleQuadrangulation.h
    LINK
      barycentricSubdivision
      dijkstra
      geometry
      triangulation
      )
else()
  set(VTK_MODULE_ENABLE_ttkMorseSmaleQuadrangulation "NO" CACHE BOOL "morseSmaleQuadrangulation module" FORCE)
  message(STATUS "morseSmaleQuadrangulation base code and wrapped disabled with new cell array layout.")
endif()

