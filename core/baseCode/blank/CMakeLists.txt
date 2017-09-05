ttk_add_base_library(blank Blank.cpp HEADERS Blank.h)
target_link_libraries(blank triangulation)
add_library(ttk::base::Blank ALIAS blank)

