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
  OPTIONAL_DEPENDS
    zfp::zfp
    ${ZLIB_LIBRARY}
  )

if (TTK_ENABLE_ZLIB)
  target_include_directories(topologicalCompression PUBLIC ${ZLIB_INCLUDE_DIR})
endif()

