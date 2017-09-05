ttk_add_base_library(bottleneckDistance BottleneckDistance.cpp Munkres.cpp
	HEADERS BottleneckDistance.h Munkres.h DataTypes.h)
target_link_libraries(bottleneckDistance triangulation persistenceDiagram)

