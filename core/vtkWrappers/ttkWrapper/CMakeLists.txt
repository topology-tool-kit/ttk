ttk_add_vtk_library(ttkWrapper ttkWrapper.cpp HEADERS ttkWrapper.h)
target_link_libraries(ttkWrapper common ttkTriangulation)

