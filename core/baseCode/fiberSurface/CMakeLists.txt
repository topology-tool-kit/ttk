option(TTK_BASE_FIBER_SURFACE_WITH_RANGE_DRIVEN_OCTREE
	"Use Range Driven Octree in Fiber Surface" ON)

# TODO: This option is kind of for backward compatibilty the way we set the
# add definition.
if (TTK_BASE_FIBER_SURFACE_WITH_RANGE_DRIVEN_OCTREE)
	add_definitions(-DwithrangeDrivenOctree)
endif()

set(DEPS geometry triangulation)

if (TTK_BASE_FIBER_SURFACE_WITH_RANGE_DRIVEN_OCTREE)
	set(DEPS ${DEPS} rangeDrivenOctree)
endif()

ttk_add_base_library(fiberSurface
	SOURCES FiberSurface.cpp
	HEADERS FiberSurface.h
	LINK ${DEPS})

