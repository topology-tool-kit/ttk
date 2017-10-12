ttk_add_base_library(contourForests
	SOURCES ContourForests.cpp
	HEADERS ContourForests.h ContourForestsTemplate.h
	LINK contourForestsTree)

option(TTK_CONTOUR_FORESTS_ENABLE_PARALLEL_SIMPLIFY
	"Enable parallel simplify for contour forests" OFF)
if (TTK_CONTOUR_FORESTS_ENABLE_PARALLEL_SIMPLIFY)
	target_compile_definitions(contourForests PUBLIC TTK_CONTOUR_FORESTS_ENABLE_PARALLEL_SIMPLIFY)
endif()

