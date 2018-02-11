find_package(Boost REQUIRED)

ttk_add_base_library(ftrGraph
   SOURCES Graph.cpp
   HEADERS FTRGraph.h FTRGraph_Template.h DynamicGraph.h Propagation.h Scalars.h 
   LINK triangulation)

target_include_directories(ftrGraph PUBLIC ${Boost_INCLUDE_DIR})

