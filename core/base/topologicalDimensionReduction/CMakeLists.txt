ttk_add_base_library(topologicalDimensionReduction 
  SOURCES
    DimensionReductionModel.cpp 
    TopologicalLoss.cpp 
    TopologicalDimensionReduction.cpp 
  HEADERS 
    DimensionReductionModel.h 
    TopologicalLoss.h 
    TopologicalDimensionReduction.h 
  DEPENDS
    persistenceDiagramWarmRestartAuction 
    ripsPersistenceDiagram)

if(TTK_ENABLE_TORCH)
target_link_libraries(topologicalDimensionReduction PUBLIC "${TORCH_LIBRARIES}")
target_include_directories(topologicalDimensionReduction PUBLIC ${TORCH_INCLUDE_DIRS})
target_compile_definitions(topologicalDimensionReduction PUBLIC TTK_ENABLE_TORCH) 
endif()
