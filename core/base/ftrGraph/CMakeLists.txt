find_package(Boost REQUIRED)

ttk_add_base_library(ftrGraph
   SOURCES Graph.cpp Mesh.cpp
   HEADERS AtomicVector.h FTRGraph.h FTRGraph_Template.h FTRGraphBFS_Template.h
           FTRGraphPrint_Template.h FTRGraphPrivate_Template.h DynamicGraph.h
           Mesh.h Propagation.h Scalars.h
   LINK triangulation unionFind)

target_include_directories(ftrGraph PUBLIC ${Boost_INCLUDE_DIR})

