if(APPLE)
  set(TTK_LTS_MACOS_DEPENDENCY "Boost::boost")
endif()

ttk_add_base_library(localizedTopologicalSimplification
  SOURCES
    LocalizedTopologicalSimplification.cpp
  HEADERS
    LocalizedTopologicalSimplification.h
    SSCPropagation.h
  DEPENDS
    triangulation
    ${TTK_LTS_MACOS_DEPENDENCY}
)
