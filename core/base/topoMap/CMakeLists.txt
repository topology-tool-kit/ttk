ttk_add_base_library(topoMap
  SOURCES
    TopoMap.cpp
  HEADERS
    TopoMap.h
  DEPENDS
    geometry
    unionFind
)

if(TTK_ENABLE_QHULL)
  if (Qhull_FOUND)
    target_link_libraries(topoMap PUBLIC Qhull::qhullcpp)
    target_link_libraries(topoMap PUBLIC Qhull::qhull_r)
    target_compile_definitions(topoMap PUBLIC Qhull_FOUND)
  endif()
  target_compile_definitions(topoMap PUBLIC TTK_ENABLE_QHULL)
endif()

if(NOT TTK_ENABLE_QHULL OR NOT Qhull_FOUND)
  target_link_libraries(topoMap PUBLIC Boost::boost)
endif()
