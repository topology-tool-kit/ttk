ttk_add_base_library(persistenceDiagramDictionary
  SOURCES
    PersistenceDiagramDictionary.cpp
    PersistenceDiagramDictionaryBorderInitializer.cpp
    PersistenceDiagramConstrainedOptimization.cpp
    InitDictRandomly.cpp
  HEADERS
    PersistenceDiagramDictionary.h
    PersistenceDiagramDictionaryBorderInitializer.h
    PersistenceDiagramConstrainedOptimization.h
    InitDictRandomly.h
  DEPENDS
    common
    persistenceDiagramAuction
    persistenceDiagramDistanceMatrix
    persistenceDiagramClustering
    
)

if(TTK_ENABLE_EIGEN)
  target_compile_definitions(persistenceDiagramDictionary PRIVATE TTK_ENABLE_EIGEN)
  target_link_libraries(persistenceDiagramDictionary PRIVATE Eigen3::Eigen)
endif()