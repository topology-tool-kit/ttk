ttk_add_base_library(topologicalOptimization
  SOURCES
    TopologicalOptimization.cpp
  HEADERS
    TopologicalOptimization.h 
  DEPENDS
    triangulation
    persistenceDiagram
    persistenceDiagramClustering
    )

if(TTK_ENABLE_TORCH)
  target_include_directories(topologicalOptimization PUBLIC ${TORCH_INCLUDE_DIRS})
  target_compile_options(topologicalOptimization PUBLIC "${TORCH_CXX_FLAGS}")
  target_link_libraries(topologicalOptimization PUBLIC "${TORCH_LIBRARIES}")
  target_compile_definitions(topologicalOptimization PUBLIC TTK_ENABLE_TORCH)
endif()
