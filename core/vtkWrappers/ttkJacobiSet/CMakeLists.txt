ttk_add_vtk_library(ttkJacobiSet ttkJacobiSet.cpp HEADERS ttkJacobiSet.h)
target_link_libraries(ttkJacobiSet jacobiSet ttkTriangulation)

