ttk_add_base_library(persistenceDiagramClustering
    SOURCES PersistenceDiagramClustering.cpp PersistenceDiagramBarycenter.cpp
    HEADERS PersistenceDiagramClustering.h PersistenceDiagramBarycenter.h
	LINK common auction persistenceDiagram bottleneckDistance kdTree)

# if (NOT (CMAKE_BUILD_TYPE MATCHES Debug))
#     set (CMAKE_BUILD_TYPE Debug)
# endif()
