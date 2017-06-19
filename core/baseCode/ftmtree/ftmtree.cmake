# if the package is not a template, uncomment the following line
ttk_add_basecode_package(triangulation)

# Boost dependencies
ttk_add_external_package(boost boost ${TTK_DIR}/third_party/boost/)

set(CMAKE_CXX_FLAGS "-D_GLIBCXX_PARALLEL")

ttk_wrapup_library(libFTMTree "FTMTree.cpp FTMTree_CT.cpp FTMTree_MT.cpp Segmentation.cpp")
