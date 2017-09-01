add_library(vtkBlank vtkBlank.cpp)
target_include_directories(vtkBlank PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:include/ttk/vtk>)
target_link_libraries(vtkBlank ttk::base::blank vtkTTKWrapper ${VTK_LIBRARIES})
add_library(ttk::vtk::blank ALIAS vtkBlank)


