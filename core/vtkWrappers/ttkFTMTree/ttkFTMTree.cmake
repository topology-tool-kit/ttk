ttk_add_baseCode_package(geometry)
ttk_add_baseCode_package(ftmtree)
ttk_add_vtkWrapper_package(ttkIdentifiers)

set(CMAKE_CXX_FLAGS "-D_GLIBCXX_PARALLEL")

ttk_add_source("ttkFTMTree.cpp")
