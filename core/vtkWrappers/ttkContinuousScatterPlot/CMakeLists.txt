ttk_add_vtk_library(ttkContinuousScatterPlot
	SOURCES ttkContinuousScatterPlot.cpp
	HEADERS ttkContinuousScatterPlot.h
	PLUGIN_VERSION 1.0
	PLUGIN_XML ContinuousScatterPlot.xml
	LINK continuousScatterPlot ttkTriangulation)
