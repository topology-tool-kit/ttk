ttk_add_base_library(persistenceDiagramDistanceMatrix
    SOURCES PersistenceDiagramDistanceMatrix.cpp PersistenceDiagramBarycenter.cpp
    HEADERS PersistenceDiagramDistanceMatrix.h PersistenceDiagramBarycenter.h
	LINK common auction persistenceDiagram bottleneckDistance kdTree)

# if (NOT (CMAKE_BUILD_TYPE MATCHES Debug))
#     set (CMAKE_BUILD_TYPE Debug)
# endif()
