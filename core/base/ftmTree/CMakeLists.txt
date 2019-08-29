ttk_add_base_library(ftmTree
  SOURCES
    FTMTree.cpp
    FTMTree_CT.cpp
    FTMTree_MT.cpp
    FTMSegmentation.cpp
  HEADERS
    FTMAtomicUF.h
      FTMAtomicVector.h
      FTMDataTypes.h
      FTMTree.h
      FTMTree_CT.h
      FTMTree_MT.h
      FTMTree_MT_Template.h
      FTMTree_Template.h
      FTMNode.h
      FTMSegmentation.h
      FTMStructures.h
      FTMSuperArc.h
  LINK
    triangulation
    geometry
    )

option(TTK_ENABLE_FTM_TREE_PROCESS_SPEED "Enable FTM tree process speed" OFF)
mark_as_advanced(TTK_ENABLE_FTM_TREE_PROCESS_SPEED)

option(TTK_ENABLE_FTM_TREE_STATS_TIME    "Enable FTM tree stats on tasks processing" OFF)
mark_as_advanced(TTK_ENABLE_FTM_TREE_STATS_TIME)

if (TTK_ENABLE_FTM_TREE_PROCESS_SPEED)
  target_compile_definitions(ftmTree PUBLIC TTK_ENABLE_FTM_TREE_PROCESS_SPEED)
endif()

# target_compile_definitions(ftmTree PUBLIC TTK_ENABLE_FTM_TREE_STATS_TIME)

