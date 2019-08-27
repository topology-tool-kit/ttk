ttk_add_base_library(persistenceDiagramsClustering
    SOURCES PersistenceDiagramsClustering.cpp PersistenceDiagramsBarycenter.cpp
    HEADERS PersistenceDiagramsClustering.h PersistenceDiagramsBarycenter.h
	LINK common auction persistenceDiagram bottleneckDistance kdTree)

# if (NOT (CMAKE_BUILD_TYPE MATCHES Debug))
#     set (CMAKE_BUILD_TYPE Debug)
# endif()
