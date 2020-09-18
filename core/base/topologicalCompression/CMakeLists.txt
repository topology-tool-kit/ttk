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
  )

if(TTK_ENABLE_ZLIB)
  target_compile_definitions(topologicalCompression PUBLIC TTK_ENABLE_ZLIB)
  target_include_directories(topologicalCompression SYSTEM PRIVATE ${ZLIB_INCLUDE_DIR})
  target_link_libraries(topologicalCompression PRIVATE ${ZLIB_LIBRARIES})
endif()

if(TTK_ENABLE_ZFP)
  target_compile_definitions(topologicalCompression PUBLIC TTK_ENABLE_ZFP)
  target_include_directories(topologicalCompression SYSTEM PRIVATE ${ZFP_INCLUDE_DIR})
  target_link_libraries(topologicalCompression PRIVATE ${ZFP_LIBRARIES})
endif()
