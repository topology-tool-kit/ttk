ttk_add_base_library(bottleneckDistance
  SOURCES BottleneckDistance.cpp Munkres.cpp
  HEADERS BottleneckDistance.h BottleneckDistanceImpl.h BottleneckDistanceMainImpl.h Munkres.h
	LINK triangulation persistenceDiagram)
