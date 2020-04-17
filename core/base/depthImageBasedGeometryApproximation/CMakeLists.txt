if ("${TTK_CELL_ARRAY_LAYOUT}" STREQUAL "SingleArray")
  ttk_add_base_library(depthImageBasedGeometryApproximation
    SOURCES
      DepthImageBasedGeometryApproximation.cpp
    HEADERS
      DepthImageBasedGeometryApproximation.h
    LINK
      triangulation
      )
else()
  set(VTK_MODULE_ENABLE_ttkDepthImageBasedGeometryApproximation "NO" CACHE BOOL "Geometry approximation module" FORCE)
  message(STATUS "Geometry approximation base code and wrapped disabled with new cell array layout.")
endif()
