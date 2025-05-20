ttk_add_base_library(topologicallyConstrainedDimensionReduction
  SOURCES
    DimensionReductionModel.cpp
    TopologicalLoss.cpp
    TopologicallyConstrainedDimensionReduction.cpp
  HEADERS
    DimensionReductionModel.h
    TopologicalLoss.h
    TopologicallyConstrainedDimensionReduction.h
  DEPENDS
    persistenceDiagramWarmRestartAuction
    ripsPersistenceDiagram
)

if (TTK_ENABLE_TORCH)
  target_link_libraries(topologicallyConstrainedDimensionReduction PUBLIC "${TORCH_LIBRARIES}")
  target_include_directories(topologicallyConstrainedDimensionReduction PUBLIC ${TORCH_INCLUDE_DIRS})
  target_compile_definitions(topologicallyConstrainedDimensionReduction PUBLIC TTK_ENABLE_TORCH)
endif()