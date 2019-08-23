option(TTK_BUILD_CONTINUOUSSCATTERPLOT_FILTER "Build the ContinuousScatterPlot filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_CONTINUOUSSCATTERPLOT_FILTER)

if(${TTK_BUILD_CONTINUOUSSCATTERPLOT_FILTER})
  ttk_register_pv_filter(ttkContinuousScatterPlot ContinuousScatterPlot.xml)
endif()
