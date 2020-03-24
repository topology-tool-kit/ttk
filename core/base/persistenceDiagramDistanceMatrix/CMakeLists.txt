ttk_add_base_library(persistenceDiagramDistanceMatrix
  SOURCES PersistenceDiagramDistanceMatrix.cpp PDDistMat.cpp
  HEADERS PersistenceDiagramDistanceMatrix.h
  LINK common auction persistenceDiagram bottleneckDistance kdTree
  )
