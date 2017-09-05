ttk_add_vtk_library(ttkFiber ttkFiber.cpp HEADERS ttkFiber.h)
target_link_libraries(ttkFiber PUBLIC common ttkTriangulation)

