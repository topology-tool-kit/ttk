if ("${TTK_CELL_ARRAY_LAYOUT}" STREQUAL "SingleArray")
  ttk_add_base_library(reebSpace
    SOURCES
      ReebSpace.cpp
    HEADERS
      ReebSpace.h
    DEPENDS
      fiberSurface
      geometry
      jacobiSet
      triangulation
      )
else()
  set(VTK_MODULE_ENABLE_ttkReebSpace "NO" CACHE BOOL "Reeb space module" FORCE)
  message(STATUS "Reeb space base code and wrapped disabled with new cell array layout.")
endif()
