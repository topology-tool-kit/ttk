ttk_add_base_library(persistenceDiagramsClustering
    SOURCES PersistenceDiagramsClustering.cpp 
    HEADERS PersistenceDiagramsClustering.h 
	LINK common persistenceDiagram bottleneckDistance kdTree persistenceDiagramsBarycenter)

# target_compile_options(persistenceDiagramsClustering PRIVATE -Wextra -pedantic)
