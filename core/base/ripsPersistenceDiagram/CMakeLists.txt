ttk_add_base_library(ripsPersistenceDiagram
  SOURCES
    RipsPersistenceDiagram.cpp
    RipsPersistenceDiagramUtils.cpp
    FastRipsPersistenceDiagram2.cpp
    PairCells.cpp
    PairCellsWithOracle.cpp
    ripserpy.cpp
  HEADERS
    RipsPersistenceDiagram.h
    RipsPersistenceDiagramUtils.h
    FastRipsPersistenceDiagram2.h
    PairCells.h
    PairCellsWithOracle.h
    ripser.h
  DEPENDS
    common
)

option(TTK_ENABLE_RIPSER_128BITS_IDS "Use 128 bits indices in Ripser" OFF)
mark_as_advanced(TTK_ENABLE_RIPSER_128BITS_IDS)

if (TTK_ENABLE_RIPSER_128BITS_IDS)
  target_compile_definitions(ripsPersistenceDiagram PUBLIC TTK_ENABLE_RIPSER_128BITS_IDS)
endif()

if(TTK_ENABLE_CGAL)
  target_compile_definitions(ripsPersistenceDiagram PUBLIC TTK_ENABLE_CGAL)
  target_link_libraries(ripsPersistenceDiagram PRIVATE ${CGAL_LIBRARIES})
endif()