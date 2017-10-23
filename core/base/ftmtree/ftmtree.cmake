# if the package is not a template, uncomment the following line
ttk_add_basecode_package(triangulation)

# Boost dependencies
ttk_add_external_header_package(boost)

ttk_wrapup_library(libFTMTree "FTMTree.cpp FTMTree_CT.cpp FTMTree_MT.cpp Segmentation.cpp")
