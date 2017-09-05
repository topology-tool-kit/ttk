add_library(blank SHARED Blank.cpp)
target_include_directories(blank PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:include/ttk/base>)
target_link_libraries(blank triangulation)
add_library(ttk::base::Blank ALIAS blank)

