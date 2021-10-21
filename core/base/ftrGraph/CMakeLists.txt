# compilation options

option(ENABLE_PROFILING "Enable GPerfTools profiling" OFF)
mark_as_advanced(ENABLE_PROFILING)

set(profiler_lib "")
if(ENABLE_PROFILING)
  set(profiler_lib profiler)
endif()

# option (TTK_DISABLE_FTR_LAZY "Disable the lazyness optimization" OFF)
# mark_as_advanced(TTK_DISABLE_FTR_LAZY)

option (TTK_ENABLE_FTR_TASK_STATS "print propagations death time (min only)" OFF)
mark_as_advanced(TTK_ENABLE_FTR_TASK_STATS)

option (TTK_ENABLE_FTR_VERT_STATS "print vertices touch" OFF)
mark_as_advanced(TTK_ENABLE_FTR_VERT_STATS)

# create lib

ttk_add_base_library(ftrGraph
  SOURCES
    Graph.cpp
    FTRSegmentation.cpp
  HEADERS
    FTRAtomicUF.h
    FTRAtomicVector.h
    FTRCommon.h
    FTRDataTypes.h
    FTRGraph.h
    FTRGraph_Template.h
    FTRGraphPrint_Template.h
    FTRGraphPrivate_Template.h
    Graph.h
    Graph_Template.h
    DynamicGraph.h
    DynamicGraph_Template.h
    FTRLazy.h
    FTRNode.h
    Mesh.h
    FTRPropagation.h
    FTRPropagations.h
    FTRScalars.h
    FTRSegmentation.h
    FTRSuperArc.h
    FTRTasks.h
  DEPENDS
    triangulation
    scalarFieldCriticalPoints
    Boost::boost
    ${profiler_lib}
    )

# target definitions

if(ENABLE_PROFILING)
  target_compile_definitions(ftrGraph PUBLIC GPROFILE=1)
endif()

# if(TTK_DISABLE_FTR_LAZY)
#    target_compile_definitions(ftrGraph PUBLIC TTK_DISABLE_FTR_LAZY)
# endif()

if(TTK_ENABLE_FTR_TASK_STATS)
  target_compile_definitions(ftrGraph PUBLIC TTK_ENABLE_FTR_TASK_STATS)
endif()

if(TTK_ENABLE_FTR_VERT_STATS)
  target_compile_definitions(ftrGraph PUBLIC TTK_ENABLE_FTR_VERT_STATS)
endif()
