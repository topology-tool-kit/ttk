ttk_add_baseCode_package(geometry)
ttk_add_baseCode_package(taskedTree)

set(CMAKE_CXX_FLAGS "-D_GLIBCXX_PARALLEL")

ttk_add_source("vtkTaskedTree.cpp")
