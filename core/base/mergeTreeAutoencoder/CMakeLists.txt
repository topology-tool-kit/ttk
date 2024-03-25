ttk_add_base_library(mergeTreeAutoencoder
  SOURCES
    MergeTreeAutoencoder.cpp
    MergeTreeAutoencoderUtils.cpp
    MergeTreeTorchUtils.cpp
  HEADERS
    MergeTreeAutoencoder.h
    MergeTreeAutoencoderUtils.h
    MergeTreeTorchUtils.h
  DEPENDS
    mergeTreePrincipalGeodesics
    geometry
)

if(TTK_ENABLE_TORCH)
  target_include_directories(mergeTreeAutoencoder PUBLIC ${TORCH_INCLUDE_DIRS})
  target_compile_options(mergeTreeAutoencoder PUBLIC "${TORCH_CXX_FLAGS}")
  target_link_libraries(mergeTreeAutoencoder PUBLIC "${TORCH_LIBRARIES}")
  target_compile_definitions(mergeTreeAutoencoder PUBLIC TTK_ENABLE_TORCH)
endif()
