ttk_add_vtk_library(ttkTriangulation ttkTriangulation.cpp ttkTriangulation.cpp ttkWrapper.cpp
	HEADERS ttkTriangulation.h ttkTriangulation.h ttkWrapper.h)
target_link_libraries(ttkTriangulation PUBLIC triangulation)

