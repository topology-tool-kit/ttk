# if the package is not a template, uncomment the following line
ttk_add_basecode_package(triangulation)

#add_library(any.hpp Boost ${TTK_DIR}/baseCode/taskedTree/lib/boost_1_64_0/boost/)
ttk_add_external_package(boost boost ${TTK_DIR}/baseCode/taskedTree/lib/boost_1_64_0/)

set(CMAKE_CXX_FLAGS "-D_GLIBCXX_PARALLEL")

ttk_wrapup_library(libTaskedTree "TaskedTree.cpp ContourTree.cpp MergeTree.cpp Segmentation.cpp")
