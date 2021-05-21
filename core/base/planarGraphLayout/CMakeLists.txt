ttk_add_base_library(planarGraphLayout
  SOURCES
    PlanarGraphLayout.cpp
  HEADERS
    PlanarGraphLayout.h
  DEPENDS
    common
)

if (TTK_ENABLE_GRAPHVIZ AND GRAPHVIZ_FOUND)
  target_compile_definitions(planarGraphLayout PRIVATE TTK_ENABLE_GRAPHVIZ)
  target_include_directories(planarGraphLayout PRIVATE ${GRAPHVIZ_INCLUDE_DIR})
  target_link_libraries(planarGraphLayout PRIVATE ${GRAPHVIZ_CDT_LIBRARY})
  target_link_libraries(planarGraphLayout PRIVATE ${GRAPHVIZ_GVC_LIBRARY})
  target_link_libraries(planarGraphLayout PRIVATE ${GRAPHVIZ_CGRAPH_LIBRARY})
  target_link_libraries(planarGraphLayout PRIVATE ${GRAPHVIZ_PATHPLAN_LIBRARY})
endif()
