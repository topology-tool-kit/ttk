ttk_add_base_library(contourForestsTree
	ContourForestsTree.cpp
	MergeTree.cpp
	Segmentation.cpp

	HEADERS
	ContourForestsTree.h
	DataTypes.h
	ExtendedUF.h
	MergeTree.h
	MergeTreeTemplate.h
	Node.h
	Segmentation.h
	Structures.h
	SuperArc.h)

target_link_libraries(contourForestsTree triangulation geometry)

