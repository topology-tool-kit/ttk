ttk_add_base_library(ripsPersistenceDiagram
  SOURCES
    RipsPersistenceDiagram.cpp
    ripserpy.cpp
  HEADERS
    RipsPersistenceDiagram.h
    ripser.h
  DEPENDS
    common
)

option(TTK_ENABLE_RIPSER_128BITS_IDS "Use 128 bits indices in Ripser" OFF)
mark_as_advanced(TTK_ENABLE_RIPSER_128BITS_IDS)

if (TTK_ENABLE_RIPSER_128BITS_IDS)
    target_compile_definitions(ripsPersistenceDiagram PUBLIC TTK_ENABLE_RIPSER_128BITS_IDS)
endif()
