find_package(Boost REQUIRED)

ttk_add_base_library(ftmTree
	SOURCES FTMTree.cpp FTMTree_CT.cpp FTMTree_MT.cpp Segmentation.cpp
	HEADERS AtomicUF.h AtomicVector.h DataTypes.h FTMTree.h FTMTree_CT.h
		FTMTree_MT.h FTMTree_MT_Template.h FTMTree_Template.h Node.h
		Segmentation.h Structures.h SuperArc.h
	LINK triangulation geometry)
target_include_directories(ftmTree PUBLIC ${Boost_INCLUDE_DIR})

option(TTK_ENABLE_FTM_TREE_PROCESS_SPEED "Enable FTM tree process speed" OFF)

if (TTK_ENABLE_FTM_TREE_PROCESS_SPEED)
	target_compile_definitions(ftmTree PUBLIC TTK_ENABLE_FTM_TREE_PROCESS_SPEED)
endif()

# target_compile_definitions(ftmTree PUBLIC TTK_ENABLE_FTM_TREE_STATS_TIME)

