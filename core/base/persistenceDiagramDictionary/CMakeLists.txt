ttk_add_base_library(persistenceDiagramDictionary
  SOURCES
    PersistenceDiagramDictionary.cpp
  HEADERS
    PersistenceDiagramDictionary.h
  DEPENDS
    common
    persistenceDiagramAuction
    constrainedGradientDescent
    persistenceDiagramClustering
    initDictPersistenceDiagram
)

if(TTK_ENABLE_EIGEN)
  target_compile_definitions(persistenceDiagramDictionary PRIVATE TTK_ENABLE_EIGEN)
  target_link_libraries(persistenceDiagramDictionary PRIVATE Eigen3::Eigen)
endif()