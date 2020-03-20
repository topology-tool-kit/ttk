ttk_add_base_library(persistenceDiagramDistanceMatrix
    SOURCES PersistenceDiagramDistanceMatrix.cpp
    HEADERS PersistenceDiagramDistanceMatrix.h
	LINK common auction persistenceDiagram bottleneckDistance kdTree)

# if (NOT (CMAKE_BUILD_TYPE MATCHES Debug))
#     set (CMAKE_BUILD_TYPE Debug)
# endif()
