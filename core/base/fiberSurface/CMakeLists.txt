option(TTK_ENABLE_FIBER_SURFACE_WITH_RANGE_OCTREE "Use Range Driven Octree in Fiber Surface" ON)
mark_as_advanced(TTK_ENABLE_FIBER_SURFACE_WITH_RANGE_OCTREE)

set(DEPS geometry triangulation)

if (TTK_ENABLE_FIBER_SURFACE_WITH_RANGE_OCTREE)
  set(DEPS ${DEPS} rangeDrivenOctree)
endif()

ttk_add_base_library(fiberSurface
  SOURCES
    FiberSurface.cpp
  HEADERS
    FiberSurface.h
  DEPENDS
    ${DEPS}
    )

if (TTK_ENABLE_FIBER_SURFACE_WITH_RANGE_OCTREE)
  target_compile_definitions(fiberSurface PUBLIC TTK_ENABLE_FIBER_SURFACE_WITH_RANGE_OCTREE)
endif()
