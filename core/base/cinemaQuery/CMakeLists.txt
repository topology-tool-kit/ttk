ttk_add_base_library(cinemaQuery
  SOURCES
    CinemaQuery.cpp
  HEADERS
    CinemaQuery.h
  DEPENDS
    triangulation
)

if(TTK_ENABLE_SQLITE3)
  target_compile_definitions(cinemaQuery PRIVATE TTK_ENABLE_SQLITE3)
  target_link_libraries(cinemaQuery PRIVATE SQLite::SQLite3)
endif()
