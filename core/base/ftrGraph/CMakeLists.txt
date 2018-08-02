ttk_add_base_library(ftrGraph
   SOURCES Graph.cpp Mesh.cpp Segmentation.cpp
   HEADERS AtomicVector.h FTRGraph.h FTRGraph_Template.h FTRGraphBFS_Template.h
           FTRGraphPrint_Template.h FTRGraphPrivate_Template.h DynamicGraph.h
           Lazy.h Mesh.h Propagation.h Propagations.h Scalars.h Segmentation.h
           LINK triangulation)

target_include_directories(ftrGraph PUBLIC ${Boost_INCLUDE_DIR})

if(OPENMP_FOUND)
   option(TTK_ENABLE_FTR_PRIORITY "Gives tasks priority, high perf improvement" OFF)

   if (OpenMP_CXX_VERSION_MAJOR GREATER_EQUAL 4 AND OpenMP_CXX_VERSION_MINOR GREATER_EQUAL 5)
      set(TTK_ENABLE_FTR_PRIORITY ON CACHE BOOL "Enable priorities on opnemp tasks" FORCE)
   endif()

   mark_as_advanced(TTK_ENABLE_FTR_PRIORITY)

   if (TTK_ENABLE_FTR_PRIORITY)
      target_compile_definitions(ftrGraph PUBLIC TTK_ENABLE_FTR_PRIORITY)
   endif()
endif()

option (TTK_ENABLE_FTR_STATS "print propagations death time" OFF)
mark_as_advanced(TTK_ENABLE_FTR_STATS)
if(TTK_ENABLE_FTR_STATS)
   target_compile_definitions(ftrGraph PUBLIC TTK_ENABLE_FTR_STATS)
endif()

option (TTK_DISABLE_FTR_LAZY "Disable the lazyness optimization" OFF)
mark_as_advanced(TTK_DISABLE_FTR_LAZY)
if(TTK_DISABLE_FTR_LAZY)
   target_compile_definitions(ftrGraph PUBLIC TTK_DISABLE_FTR_LAZY)
endif()
