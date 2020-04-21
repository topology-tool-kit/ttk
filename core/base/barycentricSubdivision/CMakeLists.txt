# if ("${TTK_CELL_ARRAY_LAYOUT}" STREQUAL "SingleArray")
  ttk_add_base_library(barycentricSubdivision
    SOURCES
      BarycentricSubdivision.cpp
    HEADERS
      BarycentricSubdivision.h
    DEPENDS
      triangulation
      geometry
    )
# else()
#   set(VTK_MODULE_ENABLE_ttkBarycentricSubdivision "NO" CACHE BOOL "barycentricSubdivision module" FORCE)
#   message(STATUS "barycentricSubdivision base code and wrapped disabled with new cell array layout.")
# endif()
