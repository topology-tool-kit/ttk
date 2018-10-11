ttk_add_base_library(bottleneckDistance
#[[<<<<<<< HEAD
  SOURCES BottleneckDistance.cpp Munkres.cpp
  HEADERS BottleneckDistance.h BottleneckDistanceImpl.h BottleneckDistanceMainImpl.h Munkres.h
	LINK triangulation persistenceDiagram kdTree)
]]
  SOURCES BottleneckDistance.cpp
  HEADERS BottleneckDistance.h BottleneckDistanceImpl.h BottleneckDistanceMainImpl.h
          Munkres.h GabowTarjan.h GabowTarjanImpl.h MatchingGraph.h
	LINK triangulation persistenceDiagram kdTree)
