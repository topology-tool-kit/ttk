ttk_add_base_library(topologicalCompression
  SOURCES
    TopologicalCompression.cpp
  HEADERS
    TopologicalCompression.h
    PersistenceDiagramCompression.h
    OtherCompression.h
  DEPENDS
    triangulation
    persistenceDiagram
    topologicalSimplification
  )

if (TTK_ENABLE_ZFP)
  target_compile_definitions(topologicalCompression PUBLIC TTK_ENABLE_ZFP)
  target_link_libraries(topologicalCompression PUBLIC zfp::zfp)
endif()
