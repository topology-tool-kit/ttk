# Allows to disable each filter
option(TTK_BUILD_IMPORTEMBEDDINGFROMTABLE_FILTER "Build the ImportEmbeddingFromTable filter" ON)
mark_as_advanced(TTK_BUILD_IMPORTEMBEDDINGFROMTABLE_FILTER)

if(${TTK_BUILD_IMPORTEMBEDDINGFROMTABLE_FILTER})
  ttk_register_pv_filter(pvImportEmbeddingFromTable ttkImportEmbeddingFromTable)
endif()
