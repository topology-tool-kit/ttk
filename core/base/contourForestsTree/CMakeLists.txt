ttk_add_base_library(contourForestsTree
	SOURCES
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
	SuperArc.h

	LINK triangulation geometry)

