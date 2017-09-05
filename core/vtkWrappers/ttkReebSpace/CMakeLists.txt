ttk_add_vtk_library(ttkReebSpace ttkReebSpace.cpp HEADERS ttkReebSpace.h)
target_link_libraries(ttkReebSpace PUBLIC reebSpace ttkTriangulation)

