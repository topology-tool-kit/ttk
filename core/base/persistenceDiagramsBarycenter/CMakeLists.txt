ttk_add_base_library(persistenceDiagramsBarycenter
	SOURCES PersistenceDiagramsBarycenter.cpp
    HEADERS PersistenceDiagramsBarycenter.h 
	LINK common persistenceDiagram bottleneckDistance kdTree)

# if (NOT (CMAKE_BUILD_TYPE MATCHES Debug))
#     set (CMAKE_BUILD_TYPE Debug)
# endif()
