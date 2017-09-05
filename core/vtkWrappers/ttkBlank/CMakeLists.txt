add_library(ttkBlank ttkBlank.cpp)
target_include_directories(ttkBlank PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:include/ttk/vtk>)
target_link_libraries(ttkBlank ttk::base::Blank ttkWrapper ${VTK_LIBRARIES})
add_library(ttk::vtk::Blank ALIAS ttkBlank)


