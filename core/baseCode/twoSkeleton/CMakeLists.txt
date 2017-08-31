ttk_add_base_library(twoSkeleton TwoSkeleton.cpp HEADERS TwoSkeleton.h)
target_link_libraries(twoSkeleton oneSkeleton zeroSkeleton)

