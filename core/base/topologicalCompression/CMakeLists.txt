ttk_add_base_library(topologicalCompression
  SOURCES
    TopologicalCompression.cpp
  HEADERS
    TopologicalCompression.h
    PersistenceDiagramCompression.h
    OtherCompression.h
  DEPENDS
    triangulation
    topologicalSimplification
    ftmTreePP
  OPTIONAL_DEPENDS
    zfp::zfp
    ZLIB::ZLIB
  )

if (TTK_ENABLE_ZLIB)
  target_include_directories(topologicalCompression PUBLIC ${ZLIB_INCLUDE_DIR})
endif()

if(TTK_ENABLE_ZFP)
  target_compile_definitions(topologicalCompression PUBLIC TTK_ENABLE_ZFP)
endif()
