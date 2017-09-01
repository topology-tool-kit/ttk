ttk_add_vtk_library(vtkTTKWrapper ttkWrapper.cpp HEADERS ttkWrapper.h)
target_link_libraries(vtkTTKWrapper common vtkTriangulation)

