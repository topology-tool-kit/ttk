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
  if(SQLite3_VERSION VERSION_GREATER_EQUAL "3.51.3")
    target_link_libraries(cinemaQuery PRIVATE SQLite3::SQLite3)
  else()
    target_link_libraries(cinemaQuery PRIVATE SQLite::SQLite3)
  endif()
endif()
