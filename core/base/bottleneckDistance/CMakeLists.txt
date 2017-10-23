ttk_add_base_library(bottleneckDistance
	SOURCES BottleneckDistance.cpp Munkres.cpp
	HEADERS BottleneckDistance.h Munkres.h DataTypes.h
	LINK triangulation persistenceDiagram)

