# if the package is not a template, uncomment the following line
 ttk_add_basecode_package(triangulation)

set(CMAKE_CXX_FLAGS "-D_GLIBCXX_PARALLEL")

 ttk_wrapup_library(libTaskedTree "TaskedTree.cpp ContourTree.cpp MergeTree.cpp Segmentation.cpp")
