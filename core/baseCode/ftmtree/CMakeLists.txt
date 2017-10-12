find_package(Boost REQUIRED)

ttk_add_base_library(ftmTree
	SOURCES FTMTree.cpp FTMTree_CT.cpp FTMTree_MT.cpp Segmentation.cpp
	HEADERS AtomicUF.h AtomicVector.h DataTypes.h FTMTree.h FTMTree_CT.h
		FTMTree_MT.h FTMTree_MT_Template.h FTMTree_Template.h Node.h
		Segmentation.h Structures.h SuperArc.h
	LINK triangulation geometry)
target_include_directories(ftmTree PUBLIC ${BOOST_INCLUDE_DIRS})

option(TTK_FTM_TREE_ENABLE_DUAL_QUEUE_COMBINE "Enable FTM tree dual queue combine" OFF)
option(TTK_FTM_TREE_ENABLE_STATS_TIME "Enable FTM tree statistics timing" OFF)
option(TTK_FTM_TREE_ENABLE_PROCESS_SPEED "Enable FTM tree process speed" OFF)

if (TTK_FTM_TREE_ENABLE_DUAL_QUEUE_COMBINE)
	target_compile_definitions(ftmTree PUBLIC TTK_FTM_TREE_ENABLE_DUAL_QUEUE_COMBINE)
endif()

if (TTK_FTM_TREE_ENABLE_STATS_TIME)
	target_compile_definitions(ftmTree PUBLIC TTK_FTM_TREE_ENABLE_STATS_TIME)
endif()

if (TTK_FTM_TREE_ENABLE_PROCESS_SPEED)
	target_compile_definitions(ftmTree PUBLIC TTK_FTM_TREE_ENABLE_PROCESS_SPEED)
endif()

