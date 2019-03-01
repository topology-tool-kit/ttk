# compilation options

option(ENABLE_PROFILING "Enable GPerfTools profiling" OFF)
mark_as_advanced(ENABLE_PROFILING)

set(profiler_lib "")
if(ENABLE_PROFILING)
  set(profiler_lib profiler)
endif()

if(OPENMP_FOUND)
   option(TTK_ENABLE_FTR_PRIORITY "Gives tasks priority, high perf improvement" OFF)

   if (OpenMP_CXX_VERSION_MAJOR GREATER_EQUAL 4 AND OpenMP_CXX_VERSION_MINOR GREATER_EQUAL 5)
      set(TTK_ENABLE_FTR_PRIORITY ON CACHE BOOL "Enable priorities on opnemp tasks" FORCE)
   endif()

   mark_as_advanced(TTK_ENABLE_FTR_PRIORITY)
endif()

# option (TTK_DISABLE_FTR_LAZY "Disable the lazyness optimization" OFF)
# mark_as_advanced(TTK_DISABLE_FTR_LAZY)

option (TTK_ENABLE_FTR_TASK_STATS "print propagations death time (min only)" OFF)
mark_as_advanced(TTK_ENABLE_FTR_TASK_STATS)

option (TTK_ENABLE_FTR_VERT_STATS "print vertices touch" OFF)
mark_as_advanced(TTK_ENABLE_FTR_VERT_STATS)

# create lib

ttk_add_base_library(ftrGraph
   SOURCES Graph.cpp Mesh.cpp Segmentation.cpp
   HEADERS AtomicVector.h FTRGraph.h FTRGraph_Template.h
           FTRGraphPrint_Template.h FTRGraphPrivate_Template.h DynamicGraph.h
           Lazy.h Mesh.h Propagation.h Propagations.h Scalars.h Segmentation.h
           LINK triangulation ${profiler_lib})

target_include_directories(ftrGraph PUBLIC ${Boost_INCLUDE_DIR})

# target definitions

if(ENABLE_PROFILING)
   target_compile_definitions(ftrGraph PUBLIC GPROFILE=1)
endif()

if (TTK_ENABLE_FTR_PRIORITY)
   target_compile_definitions(ftrGraph PUBLIC TTK_ENABLE_FTR_PRIORITY)
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
