ttk_add_base_library(skeleton
	ZeroSkeleton.cpp OneSkeleton.cpp TwoSkeleton.cpp ThreeSkeleton.cpp
	HEADERS ZeroSkeleton.h OneSkeleton.h TwoSkeleton.h ThreeSkeleton.h)
target_link_libraries(skeleton common)

